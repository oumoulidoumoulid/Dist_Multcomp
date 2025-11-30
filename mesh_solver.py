import numpy as np

class MESHSolver:
    """
    Resout le systeme d'equations MESH pour une colonne de distillation

    CONVENTION D'INDEXATION:
    ------------------------
    j = 0       : Condenseur (plateau 1)
    j = 1..N-2  : Plateaux reels
    j = N-1     : Rebouilleur (plateau N)

    FLUX (convention):
    ------------------
    L[j] : debit liquide descendant du plateau j vers j+1 (mol/s)
    V[j] : debit vapeur montant du plateau j vers j-1 (mol/s)
    x[j,:] : composition liquide au plateau j (mol fraction)
    y[j,:] : composition vapeur au plateau j (mol fraction)
    """

    def __init__(self, thermo_package, N_stages, feed_stage, P=101325):
        self.thermo = thermo_package
        self.N = int(N_stages)
        self.feed_stage = int(feed_stage)  # 1-indexed in your interface, used accordingly
        self.P = P
        self.n_comp = thermo_package.n_comp

        # State variables
        self.T = np.zeros(self.N)                    # K
        self.L = np.zeros(self.N)                    # mol/s
        self.V = np.zeros(self.N)                    # mol/s
        self.x = np.zeros((self.N, self.n_comp))     # liquid compositions
        self.y = np.zeros((self.N, self.n_comp))     # vapor compositions

        # Stream specs (set in initialize)
        self.F = None
        self.z_F = None
        self.D = None
        self.B = None
        self.R = None
        self.q = 1.0  # feed thermal condition (1 = saturated liquid)

    def initialize(self, F, z_F, D, R, q=1.0, T_feed=None):
        """
        Initialise la colonne (profils, debits, compositions).

        q : fraction liquid part of feed (q = 1 saturated liquid, q = 0 saturated vapor)
        """
        self.F = float(F)
        self.z_F = np.asarray(z_F, dtype=float)
        if self.z_F.size != self.n_comp:
            raise ValueError("z_F length must equal n_comp")
        self.D = float(D)
        self.B = self.F - self.D
        self.R = float(R)
        self.q = float(q)

        # Temperatures initiales: use feed bubble/dew for guess
        if T_feed is None:
            try:
                T_feed = self.thermo.bubble_temperature(self.P, self.z_F)
            except Exception:
                T_feed = 298.15
        try:
            T_top = self.thermo.dew_temperature(self.P, self.z_F)
        except Exception:
            T_top = T_feed - 10.0

        try:
            T_bottom = self.thermo.bubble_temperature(self.P, self.z_F)
        except Exception:
            T_bottom = T_feed + 10.0

        self.T = np.linspace(T_top, T_bottom, self.N)

        # Initialize flows using CMO separation: two constant regions
        # Above feed (rectifying): j = 0 .. feed_stage-1
        # Below feed (stripping): j = feed_stage .. N-1
        # L above = R*D ; V above = L_above + D = (R+1)*D
        # L below = R*D + (1 - q)*F ; V below = (R+1)*D + q*F
        for j in range(self.N):
            if j < (self.feed_stage):
                self.L[j] = self.R * self.D
                self.V[j] = (self.R + 1.0) * self.D
            else:
                self.L[j] = self.R * self.D + (1.0 - self.q) * self.F
                self.V[j] = (self.R + 1.0) * self.D + self.q * self.F

        # Enforce V[0] = 0 (no vapor above condenser)
        self.V[0] = 0.0

        # Initial compositions: linear profile between feed and a guessed top/bottom split
        # Use feed composition as center; slightly enrich top with lighter components (assumption)
        for i in range(self.n_comp):
            self.x[:, i] = np.linspace(min(0.99, self.z_F[i] * 1.2), max(1e-6, self.z_F[i] * 0.8), self.N)
        # Normalize x
        for j in range(self.N):
            s = np.sum(self.x[j, :])
            if s <= 0:
                self.x[j, :] = np.ones(self.n_comp) / self.n_comp
            else:
                self.x[j, :] /= s

        # compute initial y from K-values
        for j in range(self.N):
            try:
                K = self.thermo.K_values(self.T[j], self.P)
                self.y[j, :] = K * self.x[j, :]
                s = np.sum(self.y[j, :])
                if s <= 0:
                    self.y[j, :] = np.copy(self.x[j, :])
                else:
                    self.y[j, :] /= s
            except Exception:
                self.y[j, :] = np.copy(self.x[j, :])

    def solve(self, max_iter=200, tol_T=1e-3, tol_x=1e-6, relax_T=0.5, verbose=True):
        """
        Solve MESH iteratively.
        relax_T: relaxation factor for temperature update (0..1)
        """
        for it in range(1, max_iter + 1):
            T_old = self.T.copy()
            x_old = self.x.copy()
            L_old = self.L.copy()
            V_old = self.V.copy()

            # 1) K-values at current temperatures
            K = np.zeros((self.N, self.n_comp))
            for j in range(self.N):
                try:
                    K[j, :] = self.thermo.K_values(self.T[j], self.P)
                except Exception:
                    # fallback: assume ideal (K=1)
                    K[j, :] = np.ones(self.n_comp)

            # 2) Solve component balances for each component (internal nodes 1..N-2)
            # We treat x[0] and x[N-1] as boundary unknowns kept from previous iteration (Dirichlet guesses).
            # Build and solve linear system for x_component for interior stages simultaneously.
            for comp in range(self.n_comp):
                # Unknowns: x_j for j=1..N-2 (if N<=2 no interior)
                if self.N <= 2:
                    # trivial: only condenser and reboiler, compute equilibrium
                    # keep x as is and set y = K*x
                    continue

                n_int = self.N - 2
                A = np.zeros((n_int, n_int))
                b = np.zeros(n_int)

                for idx in range(n_int):
                    j = idx + 1  # actual stage index
                    # coefficients from mass balance:
                    # incoming: L[j-1]*x[j-1] + V[j+1]*y[j+1] + F_j*z_j
                    # outgoing: L[j]*x[j] + V[j]*y[j]  with y[j] = K[j]*x[j]
                    L_jm1 = self.L[j - 1]
                    L_j = self.L[j]
                    V_j = self.V[j]
                    V_jp1 = self.V[j + 1]

                    K_j = K[j, comp]
                    K_jp1 = K[j + 1, comp] if (j + 1) < self.N else K_j

                    # Build tridiagonal: a*x[j-1] + b*x[j] + c*x[j+1] = d
                    a = L_jm1
                    b_coef = -(L_j + V_j * K_j)
                    c = V_jp1 * K_jp1

                    # RHS includes feed term and known boundaries (x[0] and x[N-1])
                    F_j = self.F if (j == (self.feed_stage - 1)) else 0.0
                    z_j = self.z_F[comp] if F_j != 0 else 0.0
                    d = -F_j * z_j

                    # incorporate known x[0] when idx == 0 (j == 1) into RHS
                    if idx == 0:
                        x0 = self.x[0, comp]
                        d -= a * x0
                        a = 0.0  # moved to RHS

                    # incorporate known x[N-1] when idx == n_int-1 (j == N-2)
                    if idx == (n_int - 1):
                        xN = self.x[-1, comp]
                        d -= c * xN
                        c = 0.0

                    A[idx, idx] = b_coef
                    if idx - 1 >= 0:
                        A[idx, idx - 1] = a
                    if idx + 1 < n_int:
                        A[idx, idx + 1] = c
                    b[idx] = d

                # Solve linear system for interior x (if matrix singular, fallback to previous)
                try:
                    x_int = np.linalg.solve(A, b)
                except np.linalg.LinAlgError:
                    # fallback: small perturbation to A diagonal and retry
                    try:
                        A += np.eye(A.shape[0]) * 1e-12
                        x_int = np.linalg.solve(A, b)
                    except Exception:
                        # give up for this comp and keep old values
                        x_int = self.x[1:-1, comp].copy()

                # Store interior solution
                self.x[1:-1, comp] = x_int

            # 3) Enforce positivity and normalization on x, then compute y = K * x
            for j in range(self.N):
                # avoid negative or zero
                self.x[j, :] = np.maximum(self.x[j, :], 1e-12)
                s = np.sum(self.x[j, :])
                if s <= 0:
                    self.x[j, :] = np.ones(self.n_comp) / self.n_comp
                else:
                    self.x[j, :] /= s

            # compute y from K and x
            for j in range(self.N):
                try:
                    Krow = self.thermo.K_values(self.T[j], self.P)
                except Exception:
                    Krow = np.ones(self.n_comp)
                self.y[j, :] = Krow * self.x[j, :]
                # normalize y
                sy = np.sum(self.y[j, :])
                if sy <= 0:
                    self.y[j, :] = self.x[j, :].copy()
                else:
                    self.y[j, :] /= sy

            # 4) Update temperatures using bubble point on liquid x (with relaxation)
            for j in range(self.N):
                try:
                    T_bubble = self.thermo.bubble_temperature(self.P, self.x[j, :], self.T[j])
                    # relaxation to improve stability
                    self.T[j] = (1.0 - relax_T) * self.T[j] + relax_T * T_bubble
                except Exception as e:
                    # keep previous T and warn (no silent pass)
                    if verbose:
                        print(f"Warning: bubble temperature failed at stage {j}: {e}")

            # 5) Optionally update flows (we keep CMO assumption, but ensure V[0]=0)
            self.V[0] = 0.0

            # 6) Convergence check
            err_T = np.max(np.abs(self.T - T_old))
            err_x = np.max(np.abs(self.x - x_old))
            rel_L = np.max(np.abs((self.L - L_old) / (np.where(L_old == 0, 1e-10, L_old))))
            rel_V = np.max(np.abs((self.V - V_old) / (np.where(V_old == 0, 1e-10, V_old))))

            if verbose and (it % 10 == 0 or it == 1):
                print(f"Iter {it}: err_T={err_T:.3e}, err_x={err_x:.3e}, rel_L={rel_L:.3e}, rel_V={rel_V:.3e}")

            if err_T < tol_T and err_x < tol_x:
                if verbose:
                    print(f"Converged in {it} iterations: err_T={err_T:.3e}, err_x={err_x:.3e}")
                return True

        if verbose:
            print("Warning: solver did not converge within max_iter")
        return False

    def compute_energy_balances(self):
        """
        Compute condenser and reboiler duties.

        Returns a dict:
            condenser: {duty_W, duty_kW, duty_J_per_mol_feed, duty_kJ_per_kmol_feed, ...}
            reboiler:  { ... }
            total_duty_kW, ratio_Qr_over_Qc
        Units:
            enthalpies: J/mol
            flows: mol/s
            duties: J/s (W)
            duty per feed: J/mol ; also returns kJ/kmol numerically equal to J/mol
        """
        # Condenser (stage 0)
        T_cond = self.T[0]
        x_cond = self.x[0, :]

        if self.N > 1:
            # vapor entering condenser comes from stage 1 (below condenser)
            T_v_in = self.T[1]
            y_v_in = self.y[1, :]
            V_in = self.V[1]
        else:
            T_v_in = T_cond
            y_v_in = self.y[0, :]
            V_in = (self.R + 1.0) * self.D

        try:
            H_v_in = self.thermo.mixture_enthalpy_vapor(T_v_in, y_v_in)
        except Exception:
            H_v_in = 3.5e4  # fallback J/mol

        try:
            H_l_out = self.thermo.mixture_enthalpy_liquid(T_cond, x_cond)
        except Exception:
            H_l_out = 0.0

        L_reflux = self.L[0]
        flow_liquid_out = L_reflux + self.D
        Q_cond = flow_liquid_out * H_l_out - V_in * H_v_in  # J/s (W). typically negative

        # Reboiler (stage N-1)
        T_reb = self.T[-1]
        x_reb = self.x[-1, :]
        y_reb = self.y[-1, :]

        if self.N > 1:
            T_l_in = self.T[-2]
            x_l_in = self.x[-2, :]
            L_in_reb = self.L[-2]
        else:
            T_l_in = T_reb
            x_l_in = x_reb
            L_in_reb = self.R * self.D + self.F

        try:
            H_l_in = self.thermo.mixture_enthalpy_liquid(T_l_in, x_l_in)
        except Exception:
            H_l_in = 0.0

        try:
            H_v_out = self.thermo.mixture_enthalpy_vapor(T_reb, y_reb)
        except Exception:
            H_v_out = 3.5e4

        try:
            H_bottoms = self.thermo.mixture_enthalpy_liquid(T_reb, x_reb)
        except Exception:
            H_bottoms = 0.0

        V_out_reb = self.V[-1]
        Q_reb = V_out_reb * H_v_out + self.B * H_bottoms - L_in_reb * H_l_in  # J/s

        # per-feed metrics (protect division by zero)
        duty_cond_J_per_mol_feed = None
        duty_reb_J_per_mol_feed = None
        if self.F is not None and self.F > 0:
            duty_cond_J_per_mol_feed = Q_cond / self.F
            duty_reb_J_per_mol_feed = Q_reb / self.F

        result = {
            'condenser': {
                'duty_W': Q_cond,
                'duty_kW': Q_cond / 1000.0,
                'duty_J_per_mol_feed': duty_cond_J_per_mol_feed,
                'duty_kJ_per_kmol_feed': duty_cond_J_per_mol_feed,  # numeric equality
                'vapor_in_mol_s': V_in,
                'liquid_out_mol_s': flow_liquid_out,
                'temperature_K': T_cond,
                'temperature_C': T_cond - 273.15,
                'H_vapor_in_J_per_mol': H_v_in,
                'H_liquid_out_J_per_mol': H_l_out
            },
            'reboiler': {
                'duty_W': Q_reb,
                'duty_kW': Q_reb / 1000.0,
                'duty_J_per_mol_feed': duty_reb_J_per_mol_feed,
                'duty_kJ_per_kmol_feed': duty_reb_J_per_mol_feed,
                'liquid_in_mol_s': L_in_reb,
                'vapor_out_mol_s': V_out_reb,
                'bottoms_mol_s': self.B,
                'temperature_K': T_reb,
                'temperature_C': T_reb - 273.15,
                'H_liquid_in_J_per_mol': H_l_in,
                'H_vapor_out_J_per_mol': H_v_out,
                'H_bottoms_J_per_mol': H_bottoms
            },
            'total_duty_kW': (Q_reb + Q_cond) / 1000.0,
            'ratio_Qr_over_Qc': abs(Q_reb / Q_cond) if Q_cond != 0 else None
        }
        return result

    def get_results(self):
        energy_balance = self.compute_energy_balances()
        return {
            'stages': np.arange(1, self.N + 1),
            'temperatures': self.T,
            'compositions_liquid': self.x,
            'compositions_vapor': self.y,
            'flows_liquid': self.L,
            'flows_vapor': self.V,
            'distillate': {'composition': self.x[0, :], 'flow': self.D},
            'bottoms': {'composition': self.x[-1, :], 'flow': self.B},
            'energy_balance': energy_balance
        }
