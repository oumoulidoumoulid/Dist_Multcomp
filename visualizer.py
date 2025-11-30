import plotly.graph_objects as go
import plotly.express as px

class DistillationVisualizer:
    """
    Cree des visualisations interactives des resultats
    """

    def plot_composition_profiles(self, stages, compositions, compound_names):
        """
        Profils de composition dans la colonne
        """
        fig = go.Figure()

        colors = px.colors.qualitative.Plotly

        for i, name in enumerate(compound_names):
            # Phase liquide
            fig.add_trace(go.Scatter(
                x=compositions['liquid'][:, i],
                y=stages,
                mode='lines+markers',
                name=f'{name} (liquide)',
                line=dict(color=colors[i % len(colors)], width=2),
                marker=dict(size=6)
            ))

            # Phase vapeur
            fig.add_trace(go.Scatter(
                x=compositions['vapor'][:, i],
                y=stages,
                mode='lines+markers',
                name=f'{name} (vapeur)',
                line=dict(color=colors[i % len(colors)], width=2, dash='dash'),
                marker=dict(size=6, symbol='square')
            ))

        fig.update_layout(
            title='Profils de composition dans la colonne',
            xaxis_title='Fraction molaire',
            yaxis_title='Numero de plateau',
            yaxis=dict(autorange='reversed'),
            hovermode='closest',
            template='plotly_white'
        )

        return fig

    def plot_temperature_profile(self, stages, temperatures):
        """
        Profil de temperature
        """
        fig = go.Figure()

        fig.add_trace(go.Scatter(
            x=temperatures - 273.15,
            y=stages,
            mode='lines+markers',
            name='Temperature',
            line=dict(color='red', width=3),
            marker=dict(size=8, color='darkred')
        ))

        fig.update_layout(
            title='Profil de temperature dans la colonne',
            xaxis_title='Temperature (C)',
            yaxis_title='Numero de plateau',
            yaxis=dict(autorange='reversed'),
            template='plotly_white'
        )

        return fig
