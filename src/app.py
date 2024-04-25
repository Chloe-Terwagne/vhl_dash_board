"""
    Project :  Visualisation platform of the VHL SGE experiment
    Name : Chloé Terwagne
    Date : 20 June 2023
    Python version 3.9
"""

# IMPORT ---------------------------------------------------------------

from dash import Dash, dcc, html, Output, Input, dash_table
import dash_bootstrap_components as dbc
import plotly.express as px
import pandas as pd
from protein_3d import create_style_3d
import dash_bio as dashbio
from dash_bio.utils import PdbParser
from dash.development.base_component import Component, _explicitize_args
import plotly.graph_objs as go
import dash_daq as daq

# Launch app------------------------------------------------------------------------------------------------
app = Dash(__name__, external_stylesheets=[dbc.themes.DARKLY], suppress_callback_exceptions=True,
           meta_tags=[{'name': 'viewport', 'content': 'width=device-width, initial-scale=1.0'}])
server = app.server

# FONT & COLOR  ---------------------------------------------------------------
font_list = ["Arial", "Balto", "Courier New", "Droid Sans", "Droid Serif", "Droid Sans Mono", "Gravitas One",
             "Old Standard TT", "Open Sans", "Overpass", "PT Sans Narrow", "Raleway", "Times New Roman"]
idx_font = 0
# Color Model
light_gray = 'rgb(140, 130, 121)'
yellow = 'rgb(246, 190, 0)'
yel = "rgb(214,210,196)"
yel_exon = "rgba(214,210,196,0.1)"
dark_gray = 'rgba(41,41,41,1)'  # background
very_dark_gray = 'rgba(33,33,33,1)'  # very dark
dark_gray_transp = 'rgba(41,41,41,0.85)'
transparent = 'rgba(0,0,0,0) '

# Create color palette for conseq, clinvar, tier class and cbio
dict_cons_colors = {"Stop gained": "#e83677", "Canonical splice": "#f7a83e", "Intronic": "#8acdef",
                    "Missense": "#64bb97", "Splice site": "#243672", "Synonymous": "#4e77bb",
                    "Stop lost": "#964594",
                    "REGULATORY": "#B5A695", "5PRIME_UTR": "#B5A695", "3' UTR": "#d091bf", "REF": "#000000"}
dict_clin_colors = {"Absent ": "#979090", "Benign": "#0433FF",
                    "Conflicting interpretations of pathogenicity": "#936104",
                    "Likely benign": "#486BFC", "Likely pathogenic": "#C56B99", "Pathogenic": "#882255",
                    "Uncertain significance": "#DCA237"}
dict_cbio_color = {"Clear Cell Renal Cell Carcinoma": "#E0004B", "Papillary Renal Cell Carcinoma": "#99006B",
                   "Chromophobe Renal Cell Carcinoma": "#FFBB33", "Renal Cell Carcinoma Other": "#F580A7",
                   "Pheochromocytoma": "#00BA74", "Other": "#4F71E1", "Absent": "#D0D0D0",
                   "Renal Cell Carcinoma, NOS": "#F580A7", "Pancreatic Neuroendocrine": "#CDF800"}
dict_tier_class_green_red = {"Neutral": 'rgb(143,153,62)', "Intermediate": "#d7dc99", 'LOF2': '#f3a66e',
                             'LOF1': 'rgb(147,39,44)'}
dict_tier_class_c_blind_friendly = {"Neutral": 'rgb(143,153,62)', "Intermediate": '#e7f7d5', 'LOF2': '#d091bb',
                                    'LOF1': 'rgb(80, 7, 120)'}
CUSTOM_CAT_ORDER = ["Neutral", "Intermediate", 'LOF2', 'LOF1'] + ['Synonymous', 'Missense', "Intronic",
                                                                  "3' UTR", "Splice site", "Missense",
                                                                  "Canonical splice", "Stop lost", "Stop gained"] + [
                       "NA", "Absent ", "Benign", "Likely benign", "Uncertain significance",
                       "Conflicting interpretations of pathogenicity", "Likely pathogenic", "Pathogenic"] + ["Absent",
                                                                                                             "Other",
                                                                                                             "Clear Cell Renal Cell Carcinoma",
                                                                                                             "Papillary Renal Cell Carcinoma",
                                                                                                             "Chromophobe Renal Cell Carcinoma",
                                                                                                             "Renal Cell Carcinoma Other",
                                                                                                             "Pheochromocytoma",
                                                                                                             "Renal Cell Carcinoma, NOS",
                                                                                                             "Pancreatic Neuroendocrine"]
DICT_COL_REG = {**dict_cons_colors, **dict_clin_colors, **dict_cbio_color, **dict_tier_class_green_red}
DICT_COL_BLIND = {**dict_cons_colors, **dict_clin_colors, **dict_cbio_color, **dict_tier_class_c_blind_friendly}

# glossary padding
cell_style = {'padding-bottom': '20px', 'font-weight': 'bold', 'color': yel}  # more title type
bottom_style = {'padding-bottom': '15px'}  # more cell type

# pd display
templates = 'plotly_dark'
pd.set_option('display.width', 900)
pd.set_option('display.max_columns', 350)
pd.set_option("display.max_rows", None)


# hover display
# define hover for all
hover_columns = ['variant_id', 'cHGVS', 'pHGVS', 'consequence','function_score_final', 'tier_class']
hover_text = ["<b>%{customdata[0]}</b>",
              "cHGVS: %{customdata[1]}",
              "pHGVS: %{customdata[2]}",
              "Consequence: %{customdata[3]}",
              "SGE function score: %{customdata[4]:.2f}",
              "SGE function class: %{customdata[5]}"]


# FUNCTION & CLASS ----------------------------------------------------------------------------------------------------

# class BooleanSwitch(Component) copied from https://github.com/plotly/dash-daq/blob/master/dash_daq/BooleanSwitch.py
class BooleanSwitch(Component):
    @_explicitize_args
    def __init__(self, id=Component.UNDEFINED, on=Component.UNDEFINED, color=Component.UNDEFINED,
                 vertical=Component.UNDEFINED, disabled=Component.UNDEFINED, theme=Component.UNDEFINED,
                 label=Component.UNDEFINED, labelPosition=Component.UNDEFINED, className=Component.UNDEFINED,
                 style=Component.UNDEFINED, persistence=Component.UNDEFINED, persisted_props=Component.UNDEFINED,
                 persistence_type=Component.UNDEFINED, size=Component.UNDEFINED, **kwargs):
        self._prop_names = ['id', 'className', 'color', 'disabled', 'label', 'labelPosition', 'on', 'persisted_props',
                            'persistence', 'persistence_type', 'size', 'style', 'theme', 'vertical']
        self._type = 'BooleanSwitch'
        self._namespace = 'dash_daq'
        self._valid_wildcard_attributes = []
        self.available_properties = ['id', 'className', 'color', 'disabled', 'label', 'labelPosition', 'on',
                                     'persisted_props', 'persistence', 'persistence_type', 'size', 'style', 'theme',
                                     'vertical']
        self.available_wildcard_properties = []
        _explicit_args = kwargs.pop('_explicit_args')
        _locals = locals()
        _locals.update(kwargs)  # For wildcard attrs
        args = {k: _locals[k] for k in _explicit_args if k != 'children'}
        for k in []:
            if k not in args:
                raise TypeError(
                    'Required argument `' + k + '` was not specified.')
        super(BooleanSwitch, self).__init__(**args)


def color_bar_structure():
    # Create a scatterplot with invisible points
    fig_color_bar = px.scatter(
        df, x="rna_score", y='rna_score',
        color=df['average_fs_missense_at_aa_rna'],
        color_continuous_scale=['#DE2A17', '#823B6F', '#38378E'],  # red to blue
    )

    # Hide the points by setting opacity to 0 and marker size to 0
    fig_color_bar.update_traces(
        opacity=0,
        marker=dict(size=0)
    )

    fig_color_bar.update_layout(
        plot_bgcolor=transparent,
        paper_bgcolor=transparent,
        xaxis=dict(showgrid=False, visible=False), yaxis=dict(showgrid=False, visible=False),  # Hide axis
        height=280,
        width=180,
        font_color=yel,
        font_family=font_list[idx_font],
        coloraxis_colorbar=dict(
            title=None,
            lenmode="pixels",
            len=150,
            yanchor="top",
            y=0.5,
            x=-2,
            borderwidth=2,
            bordercolor=yel,
            tickcolor=yel,
            tickvals=[df['average_fs_missense_at_aa_rna'].max(), df['average_fs_missense_at_aa_rna'].min()],
            tickwidth=2,
            tickmode='array',
            ticks="outside",
            ticktext=["Neutral", "LoF"],
            ticklabelposition='inside left'
        )
    )
    return fig_color_bar


def variant_first_search_dropdown(list_var, df):
    print('list_var', list_var)
    # variant to have first in the dropdown
    mask = df['cHGVS'].isin(list_var)
    df_display_first = df[mask]

    # Set the 'cHGVS' column as the index
    df_display_first.set_index('cHGVS', inplace=True)

    # Reindex the DataFrame based on the desired order in list_var_to_display_first
    df_display_first = df_display_first.reindex(list_var_to_display_first)

    # Reset the index to restore the original DataFrame structure
    df_display_first.reset_index(inplace=True)
    df_remaining = df[~mask]
    list_all_var_key = [{'label': f'{row["cHGVS"]}', 'value': row['variant_id']} for index, row in
                        df_display_first.iterrows()] + \
                       [{'label': f'{row["cHGVS"]}', 'value': row['variant_id']} for index, row in
                        df_remaining.iterrows()]

    return list_all_var_key


def get_structure_file(selected_pdb_file):
    if selected_pdb_file == ['VHL_B_H_C']:
        selected_pdb_file = '1LM8_vbch_isolated.pdb?raw=true'
    else:
        selected_pdb_file = '1LM8_vhl_isolated.pdb?raw=true'
    #old dir
    #parser = PdbParser('https://github.com/Chloe-Terwagne/board_materials/blob/main/input/3d_structure/' + selected_pdb_file)
    print("previous path:", 'https://github.com/Chloe-Terwagne/board_materials/blob/main/input/3d_structure/' + selected_pdb_file)
    print("new path:", 'https://github.com/Chloe-Terwagne/vhl_dash_board/tree/main/src/assets/input/3d_structure/' + selected_pdb_file)
    parser = PdbParser('https://github.com/Chloe-Terwagne/vhl_dash_board/tree/main/src/assets/input/3d_structure/' + selected_pdb_file)    
    parser.mol3d_data()
    return parser.mol3d_data()


def print_var_score_for_selected_residue(df, aa_name):
    """
    Print all variant at the residue
    """
    if len(list(df['variant_id'])) < 1:
        text = html.Div(
            [html.Br(), html.Div(aa_name),
             html.Div("No missense variants with RNA score ≥ -2 correspond to this residue"), html.Br()])
        return text

    else:
        if len(list(df['variant_id'])) == 1:
            variant_nb = 'This residue has ' + str(
                len(list(df['variant_id']))) + ' missense variant with RNA score ≥ -2.' + '\n'
        if len(list(df['variant_id'])) > 1:
            variant_nb = 'This residue has ' + str(len(list(df[
                                                                'variant_id']))) + ' missense variants with RNA score ≥ -2 with an average function score of ' + str(
                round(list(df['average_fs_missense_at_aa_rna'])[0], 2)) + '.\n'

        variant_table = html.Table(
            [html.Tr([
                html.Th('Variant', style={'padding': '25px'}),
                html.Th('cHGVS', style={'padding': '25px'}),
                html.Th('pHGVS', style={'padding': '25px'}),
                html.Th('SGE Function Score', style={'padding': '25px'})
            ])] +
            [html.Tr([
                html.Td('Variant ' + str(i + 1), style={'padding-top': '6px', 'padding-right':'25px','padding-left':'25px' }),
                html.Td(df['cHGVS'].iloc[i], style={'padding-top': '6px', 'padding-right':'25px','padding-left':'25px' }),
                html.Td(df['pHGVS'].iloc[i], style={'padding-top': '6px', 'padding-right':'25px','padding-left':'25px' }),
                html.Td(round(df['function_score_final'].iloc[i], 2), style={'padding-top': '6px', 'padding-right':'25px','padding-left':'25px' })
            ]) for i in range(len(df))],
            style={'borderCollapse': 'collapse', 'border': '1px solid ' + light_gray})
        text = html.Div([
            html.Br(),
            html.Div(aa_name),
            html.Div(variant_nb),
            html.Br(),
            html.Br(),
            variant_table,
        ])

        return text


# MAIN ---------------------------------------------------------------------------------------------------------------
# data
df = pd.read_csv('https://github.com/Chloe-Terwagne/vhl_dash_board/blob/main/src/assets/input/vhl_preprocess_df.csv?raw=true')
df['pHGVS'] = df['pHGVS'].fillna('N/A')
df['consequence'] = df['consequence'].replace('Non synonymous', 'Missense')
print(df.consequence.value_counts())
exon_dict = {'exon 1b': [10141958, 10142087], 'exon 1a': [10142075, 10142202], 'exon 1p': [10142743, 10142876],
             'exon 2': [10146499, 10146644], 'exon 3a': [10149760, 10149887], 'exon 3b': [10149868, 10150002]}
# Get text
github_link = html.Div([
    html.A(
        id='gh-link',
        children=['View on GitHub'],
        href="https://github.com/Chloe-Terwagne/vhl_dash_board",
        style={'margin-left': '20px', 'color': 'white', 'border': 'white', "text-decoration": 'none'},
        target="_blank"
    ),
    html.Img(src='https://github.com/Chloe-Terwagne/vhl_dash_board/blob/main/src/assets/GitHub-Mark-64px.png?raw=true',

                      style={'height': '30px', 'margin-left': '10px', 'margin-right': '20px'},
             ),
], style={
    'display': 'flex',
    'align-items': 'center',
    'justify-content': 'flex-start',
    'padding': '15px',
    'background': very_dark_gray,
    'color': 'white',
    'height': '50px',
    'width': 'auto',
    'border': 'solid 0.5px white',
})

links_content = dbc.Card(
    [html.Div(
        dbc.CardImg(
            src="https://github.com/Chloe-Terwagne/vhl_dash_board/blob/main/src/assets/CRICK_Logotype_white_RGB.png?raw=true",
            style={'width': '100px'},  # Set the specific width here
            top=True
        ),
        style={'position': 'absolute', 'top': '-20px', 'right': '-20px'}  # Positioning the image at the top right
    ),

        dbc.CardBody(
            [
                html.Div([
                    html.Div([
                        github_link
                    ], style={'display': 'flex', 'flex-direction': 'row', 'justify-content': 'flex-start', 'margin-left': '300px'}),
                    html.Div([
                        dbc.CardLink(["M.Buckley", html.Em(" et al."), ', 2023'],
                                     href="https://doi.org/10.1101/2023.06.10.542698",
                                     target="_blank", className='custom-link', style={'text-decoration': 'none'}),
                    ], style={'display': 'flex', 'flex-direction': 'column', 'justify-content': 'flex-start'}),
                    html.Div([
                        dbc.CardLink("Genome Function Lab website", href="https://www.crick.ac.uk/research/labs/greg-findlay/",
                                     target="_blank", className='custom-link', style={'text-decoration': 'none'}),
                    ], style={'display': 'flex', 'flex-direction': 'column', 'justify-content': 'flex-start',
                              'margin-right': '300px'})
                ], style={'display': 'flex', 'justify-content': 'space-between', 'align-items': 'center'}),
                html.Br(),
            ]
        )
    ],
    style={"height": "100px", 'background-color': transparent, 'border': 'solid 2px',
           'border-radius': '20px',
           'border-color': light_gray,
           'overflow': 'hidden'}
)

# Build your components------------------------------------------------------------------------------------------------
# 3D parsing & styling
parser = PdbParser('https://github.com/Chloe-Terwagne/vhl_dash_board/blob/main/src/assets/input/3d_structure/1LM8_vhl_isolated.pdb?raw=true')

v_data = parser.mol3d_data()
styles = create_style_3d(
    df, 'average_fs_missense_at_aa_rna', v_data['atoms'], visualization_type='cartoon', color_element='residue_score')
vhl_3D = dashbio.Molecule3dViewer(id='dashbio-default-molecule3d', modelData=v_data, styles=styles, backgroundOpacity=0,
                                  selectionType='residue', backgroundColor="black", height=600,
                                  width=735)  # ,width=735)  # , zoom=dict(factor=1.9,animationDuration=30000, fixedPath=False))

overview_title = dcc.Markdown(children='', style=dict(font_family=font_list[idx_font], font_color=yel))
list_var_to_display_first = ['c.500G>A','c.233A>G','c.292T>C','c.473T>C',
                             'c.351G>T','c.194C>G','c.484T>C','c.334T>A','c.351G>T']

variant_highlight_dropd = dcc.Dropdown(options=variant_first_search_dropdown(list_var_to_display_first, df), multi=True,
                                       placeholder="Select or type variant(s) to highlight",
                                       className='my-custom-dropdown', style={'z-index': '2'})
var_table = dash_table.DataTable(data=[], columns=[], style_table={'overflowX': 'auto', 'backgroundColor': dark_gray},
                                 cell_selectable=False,
                                 # Background color
                                 style_data={'color': yel},  # Font color for data cells
                                 style_header={'backgroundColor': very_dark_gray, 'color': yellow},  # Header style
                                 style_cell={
                                     'backgroundColor': dark_gray_transp,  # Background color for cells
                                     'border': '1px solid white'},  # Border color
                                 )
overview_display = dcc.RadioItems(options=["SGE Function Score", "Variants expanded by nucleotide type"],
                                  value='SGE Function Score', labelClassName="custom-text p-3", labelStyle={'display': 'inline-block'},
                                  style={"margin-right": "0px!important", 'padding': '0px!important'})
overview_dropdown = dcc.Dropdown(options=[
    {'label': 'ClinVar', 'value': 'clinvar_simple'},
    {'label': 'Consequence', 'value': 'consequence'},
    {'label': 'Function Class', 'value': 'tier_class'},
    {'label': 'Cancer Type', 'value': 'Cancer_type_single'}
],
    placeholder="Select color category", value='consequence',
    clearable=False, className='my-custom-dropdown')
at_scale = BooleanSwitch(on=False, size=25, label=dict(label="Genomic position at scale", style=dict(font_color=yel)),
                         color=yellow, labelPosition="left",
                         style={"margin-right": "0px", "margin-top": "0px", "margin-bottom": "-80px", 'padding': '0px'})
overview_graph = dcc.Graph(figure={}, config={'staticPlot': False, 'scrollZoom': False, 'doubleClick': 'reset',
                                              'showTips': True, 'displayModeBar': 'hover', 'displaylogo': False,
                                              'modeBarButtonsToRemove': ['lasso2d', 'zoomIn2d', 'zoomOut2d',
                                                                         'autoScale2d']},
                           style={"margin-top": "-45px", "margin-bottom": "-75px", 'padding': '0px'}, selectedData=None)
color_blind_option = BooleanSwitch(on=False, size=25,
                                   label=dict(label="Color blind friendly", style=dict(font_color=yel)),
                                   color='rgb(80, 7, 120)', labelPosition="left")
two_d_graph = dcc.Graph(figure={},
                        config={'staticPlot': False, 'scrollZoom': False, 'doubleClick': 'reset', 'showTips': True,
                                'displayModeBar': 'hover', 'displaylogo': False,
                                'modeBarButtonsToRemove': ['lasso2d', 'select2d',
                                                           'autoScale2d'], 'watermark': False})
two_d_title = dcc.Markdown(children='all variant')
y_dropdown = dcc.Dropdown(options=[
                          {'label': 'CADD phred', 'value': 'CADD.phred'},
                          {'label': 'VARITY', 'value': 'VARITY_R'},
                          {'label': 'REVEL', 'value': 'REVEL'},
                          {'label': 'SpliceAI', 'value': 'max_spliceAI'}],
                          value='CADD.phred', clearable=False, className='my-custom-dropdown')
x_dropdown = dcc.Dropdown(options=[{'label': 'SGE Function Score', 'value': 'function_score_final'},
                                   {'label': 'RNA score', 'value': 'rna_score'}],
                          value='function_score_final', clearable=False, className='my-custom-dropdown')
mol_viewer_colorbar = dcc.Graph(figure=color_bar_structure(),  # which return fig_color_bar,
                                config={'staticPlot': True, 'scrollZoom': False, 'showTips': False,
                                        'displayModeBar': False, 'watermark': False}, style={"margin-top": '-180px'})
pdb_selector_drop = dcc.Checklist(id='pdb-selector',
                                  options=[{'label': 'Add ELOC, ELOB and HIF 1A', 'value': 'VHL_B_H_C'}],
                                  labelClassName="custom-text p-3",
                                  style={'position': 'relative', "bottom": "-103px", "margin": "0px", "padding": "0px"})
vizua_type_3d = dcc.RadioItems(id='vizua_type_3d', options={'sphere': 'Sphere', 'cartoon': 'Cartoon', 'stick': 'Stick'},
                               value='sphere', labelClassName="custom-text p-3", labelStyle={'display': 'inline-block'},
                               style={'position': 'relative', "bottom": "-50px", "margin": "0px", "padding": "0px"})

# Customize Layout--------------------------------------------------------------------------------------------
row_style = {'display': 'flex', 'flex-wrap': 'wrap', 'align-items': 'stretch'}
app.layout = \
    dbc.Container([
        dbc.Row([html.Br()]),
        dbc.Row([
            dbc.Col(html.H1("Saturation Genome Editing of VHL", className='custom-h1'), width={'size': 7, 'offset': 2}, ),
            dbc.Col([
                dbc.Row(color_blind_option, className="my-custom-switch")], width={'size': 2}, align='right')
        ], justify='between'),
        dbc.Row([html.Br()]),
        dbc.Row([
            dbc.Col([variant_highlight_dropd], width={'size': 4}),
            dbc.Col([overview_dropdown], width={'size': 2}),
        ], justify='between'),
        dbc.Row([var_table]),
        dbc.Row([html.Br()]),
        dbc.Row([html.Br()]),

        dbc.Row([dbc.Col(overview_display),
                 ], justify='between'),
        dbc.Row([
            dbc.Col(overview_graph, width=12)
        ], justify='around'),
        dbc.Row([dbc.Col([at_scale], className="my-custom-switch", width={'size': 2, 'offset': 10})]),
        # Combined Graph 2 and Graph 3 ----------------------
        dbc.Row(dbc.Col([pdb_selector_drop], width={'size': 2, 'offset': 10})),
        dbc.Row(dbc.Col([vizua_type_3d], width={'size': 3, 'offset': 6})),
        dbc.Row(
            [
                # Graph 2
                dbc.Col(
                    [
                        two_d_graph,
                        dbc.Row(
                            [
                                dbc.Col(
                                    html.Div([
                                        html.H5("Select a x-axis experimental score"),
                                        x_dropdown
                                    ], className='my-custom-title',
                                        style={'background-color': dark_gray, 'padding': '20px 20px 20px 20px',
                                               'margin': '0'}),
                                    md=6  # 6 columns for x_dropdown
                                ),
                                dbc.Col(
                                    html.Div([
                                        html.H5('Select a y-axis predictor score'),
                                        y_dropdown
                                    ], className='my-custom-title',
                                        style={'background-color': dark_gray, 'padding': '20px 20px 20px 20px',
                                               'margin-left': '-30px'}),
                                    md=6
                                ),
                            ],
                        ),
                    ],
                    md=6
                ),

                # Graph 3
                dbc.Col(
                    [
                        dbc.Row(vhl_3D),  # 3D protein
                        dbc.Row(dbc.Col([mol_viewer_colorbar], md=6)),
                        dbc.Row(dbc.Col(
                            html.H1("Averaged missense variant function score per residue mapped on VHL structure",
                                    className='custom-h1', style={
                                    'font-size': '18pt', 'text-align': 'left', "margin-left": '45px',
                                    "margin-top": '-80px', 'position': 'relative'}), width={'size': 9, 'offset': 2})),
                        dbc.Row([
                            dbc.Col(
                                [
                                    html.Div(
                                        id='default-molecule3d-output',
                                        style={
                                            'background-color': dark_gray_transp,
                                            'padding': '15px',
                                            'padding-bottom': '102px',
                                            'position': 'relative',
                                            "margin-top": '0px',
                                            "margin-left": '2px',
                                            'z-index': '1'
                                        }
                                    ),
                                ]
                            )

                        ],
                            className='custom-text_left',
                        ),
                    ],
                ),
            ],
        ),
        dbc.Row(html.Br()),
        dbc.Row(links_content),
        dbc.Row([
            dbc.Col(
                html.P(html.Small([
                    "Coded with care by ",
                    html.A("Chloé Terwagne", href="https://twitter.com/ChloeTerwagne")
                ])),
                style={'text-align': 'right'}
            )
        ]),
    ], fluid=True)


# Callback --------------------------------------------------------------------------------
@app.callback(
    Output(var_table, 'data'),
    Output(var_table, 'columns'),
    Input(variant_highlight_dropd, 'value')
)
def update_datatable(selected_variants):
    # If no variants are selected, show an empty DataTable
    if selected_variants is None or selected_variants == []:
        data, col = [], []
        return data, col

    # Filter the DataFrame based on selected variants
    subset_df = df[df['variant_id'].isin(selected_variants)]

    # Create DataTable data and columns from the subset_df
    data = subset_df.to_dict('records')
    col = [
        {'name': 'Variant', 'id': 'variant_id'},
        {'name': 'cHGVS', 'id': 'cHGVS'},
        {'name': 'pHGVS', 'id': 'pHGVS'},
        {'name': 'Consequence', 'id': 'consequence'},
        {'name': 'Function Class', 'id': 'tier_class'},
        {'name': 'Function Score', 'id': 'function_score_final', 'type': 'numeric', 'format': {'specifier': '.2f'}},
        {'name': 'RNA score', 'id': 'rna_score', 'type': 'numeric', 'format': {'specifier': '.2f'}},
        {'name': 'ClinVar', 'id': 'clinvar_simple'}
    ]
    return data, col


@app.callback(
    Output(component_id=overview_graph, component_property='figure'),
    Input(overview_dropdown, 'value'),
    Input(overview_display, 'value'),
    Input(color_blind_option, 'on'),
    Input(at_scale, 'on'),
    Input(variant_highlight_dropd, 'value')
)
def update_overview_graph(column_name, y_axis_nucleotide, color_blind, at_scale, variant_highlight):
    df_temp = df
    if at_scale:
        x_overv = 'hg38_pos'
    else:
        x_overv = 'index'

    # Get transparency if variant selected
    if variant_highlight is None or variant_highlight == []:
        transparency, mark_size, marker_line_width, ref_col = 1, 8, 0, yellow
    else:
        transparency, mark_size, marker_line_width, ref_col = 0.45, 8, 3, yel

    if color_blind:
        colors = DICT_COL_BLIND
    else:
        colors = DICT_COL_REG

    height_grph, marker_symb = 340, "circle"
    limit = (-3.745898895, 0.590402626)
    y_axis = 'function_score_final'
    yaxis_dict = dict(showgrid=True, gridcolor=yel_exon, visible=True, zeroline=False, linecolor=None, linewidth=1,
                      title='Function Score')
    xaxis_dict = dict(showgrid=False, visible=True, zeroline=False, linecolor=None, linewidth=1, showticklabels=False,
                      title="Unscaled Genomic Position")

    if y_axis_nucleotide == "Variants expanded by nucleotide type":
        marker_symb, y_axis, height_grph = "square", 'alt_pos', 340
        yaxis_dict = dict(showgrid=False, zeroline=False, title='Nucleotide',
                          tickvals=[-3.745898895, -2.3004650546666667, -0.8550312143333332, 0.590402626],
                          ticktext=['T', 'G', 'C', 'A'])

    fig = go.Figure()
    # Iterate through the unique categories in your data in the custom order
    for category in CUSTOM_CAT_ORDER:
        if category in df_temp[column_name].unique():
            category_data = df_temp[df_temp[column_name] == category]
            scatter_trace = go.Scatter(
                x=category_data[x_overv],
                y=category_data[y_axis],
                mode='markers',
                customdata=category_data[hover_columns],
                marker=dict(size=mark_size, symbol=marker_symb, color=colors[category], opacity=transparency),
                hovertemplate="<br>".join(hover_text),
                name=category)

            # Add the Scatter trace to the figure
            fig.add_trace(scatter_trace)

    # re-plot highlighted variants
    if variant_highlight is not None and variant_highlight != []:
        subset_var_highlight_df = df_temp[df_temp['variant_id'].isin(variant_highlight)]
        print("subset_var_highlight_df:")
        print(subset_var_highlight_df)
        highlight_trace = go.Scatter(
            x=subset_var_highlight_df[x_overv],
            y=subset_var_highlight_df[y_axis],
            mode='markers',
            customdata=subset_var_highlight_df[hover_columns],
            marker=dict(size=mark_size + 2, symbol=marker_symb, line=dict(width=marker_line_width, color=yellow),
                        color=[colors[key] for key in subset_var_highlight_df[column_name]]),
            hovertemplate="<br>".join(hover_text),
            name="Highlighted variant",
            opacity=1,
        )
        fig.add_trace(highlight_trace)

    # Add the scatter trace for the reference variant if applicable
    if y_axis_nucleotide == "Variants expanded by nucleotide type":
        ref_trace = go.Scatter(
            x=df_temp[x_overv],
            y=df_temp['ref_pos'],
            mode='markers',
            marker=dict(
                size=mark_size,
                symbol="square-open",
                color=ref_col,
                opacity=transparency
            ),
            hovertemplate="<br>".join(["<b>Reference allele</b>"]),
            name=''
        )
        fig.add_trace(ref_trace)
        # Add shape for intron
        if at_scale:
            start, end = 10141958 - 200, 10150002 + 200
            intron_shape = go.layout.Shape(
                type='line',
                x0=start,
                x1=end,
                y0=limit[0] + (limit[1] - limit[0]) / 2,
                y1=limit[0] + (limit[1] - limit[0]) / 2,
                line=dict(color=yel_exon, width=2),
                layer='below'
            )
            fig.add_shape(intron_shape)

    if y_axis_nucleotide == "SGE Function Score":
        # add zeroline
        fig.update_yaxes(zeroline=True, zerolinewidth=2, zerolinecolor=yel_exon, ticks="outside",
                         tickcolor='rgb(255,255,255)')

    # add exon
    if at_scale:
        xaxis_dict = dict(showgrid=False, visible=True, zeroline=False, linecolor=None, linewidth=1,
                          title="Genomic position")
        for k in exon_dict:
            start, end = exon_dict[k][0] - 0.5, exon_dict[k][1] + 0.5
            enlarge = 0.8
            texex = k
            pos_xanchor = 'center'
            fig.add_shape(type="rect",
                          x0=start, y0=limit[0] - enlarge, x1=end, y1=limit[1] + enlarge,
                          line=dict(color=yel_exon, width=2),
                          fillcolor=yel_exon, layer='below')

            fig.add_annotation(x=start, y=limit[1] + 0.1,
                               text=texex,
                               showarrow=False,
                               font=dict(color=yel, family=font_list[idx_font], size=10), textangle=270,
                               xanchor=pos_xanchor, yanchor='bottom', bgcolor=dark_gray, opacity=1)
        fig.add_shape(type="rect",
                      x0=start, y0=limit[0] - enlarge - 0.2, x1=end, y1=limit[1] + enlarge,
                      line=dict(color=transparent, width=2),
                      fillcolor=transparent, layer='below')
    dict_label_leg = {'clinvar_simple': 'ClinVar', 'consequence': 'Consequence', 'tier_class': 'Function Class', 'Cancer_type_single': 'Cancer Type'}

    fig.update_layout(
        plot_bgcolor=transparent,
        paper_bgcolor=transparent,
        xaxis=xaxis_dict,
        yaxis=yaxis_dict,
        font_family=font_list[idx_font],
        legend=dict(orientation='h', yanchor='top', y=1.2, xanchor='left', x=0,
                    title="<b>" + dict_label_leg[column_name] + ' Annotation<b>', font_family=font_list[idx_font]),
        font_color=yel,
        modebar=dict(
            bgcolor=transparent,
            activecolor=yel,
            color=yellow)
    )

    return fig.update_layout(uirevision=True)


@app.callback(
    Output(component_id=two_d_graph, component_property='figure'),
    Input(overview_dropdown, 'value'),
    Input(component_id=overview_graph, component_property="selectedData"),
    Input(x_dropdown, 'value'),
    Input(y_dropdown, 'value'),
    Input(variant_highlight_dropd, 'value'),
    Input(color_blind_option, 'on')
)
def update_2d_graph(color_column, slct_data, x_col, y_col, highlight_var, color_blind):
    black3dbg = dict(showgrid=True, gridcolor=yel_exon, gridwidth=0.5,
                     zeroline=False)

    if color_blind:
        colors = DICT_COL_BLIND
    else:
        colors = DICT_COL_REG

    df_t = df
    fig2 = go.Figure()

    #  selection with no points inside
    if slct_data is not None and slct_data['points'] == []:
        empty_trace = go.Scatter()
        fig2.add_trace(empty_trace)
        title = "Please select at least one variant"

    else:
        # if subset of point( >< not all points)
        if (slct_data is not None) and slct_data != {'points': []}:
            # Remove reference allele from slct_data -> remove points that doesn't have custom data columns
            slct_data['points'] = [point for point in slct_data['points'] if 'customdata' in point]

            if not slct_data['points']:
                empty_trace = go.Scatter()
                fig2.add_trace(empty_trace)
                title = "Please select at least one variant"
            else:
                # subset data based on selection
                var = [slct_data['points'][i]['customdata'][0] for i in range(len(slct_data['points']))]

                print("var:", var)
                df_t = df_t[df_t.variant_id.isin(var)]
                title = "Variants selected"
        else:
            subtittle = "<br><sup>Choose the rectangle tool in the menu bar of the gene overview above to subset variants of interest.</sup>"
            title = "All variants" + subtittle

        # highlighted variants settings
        if highlight_var is not None and highlight_var != []:
            transparency = 0.45
            # Create a DataFrame for highlighted points
            subset_var_highlight_df = df_t[df_t['variant_id'].isin(highlight_var)]
        else:
            transparency = 1
            subset_var_highlight_df = pd.DataFrame()

        # Iterate through the unique categories in your data in the custom order
        for category in CUSTOM_CAT_ORDER:
            # Iterate through the unique categories in your data
            if category in df_t[color_column].unique():
                # Filter data for the current category
                category_data = df_t[df_t[color_column] == category]

                # Create a Scatter trace for the current category
                scatter_trace = go.Scatter(
                    x=category_data[x_col],
                    y=category_data[y_col],
                    mode='markers',
                    customdata=category_data[hover_columns],
                    marker=dict(
                        size=6,
                        color=colors[category],  # Use color from dict_color_consq
                        opacity=transparency
                    ),
                    hovertemplate="<br>".join(hover_text),
                    name=category  # Set the name for the legend
                )

                # Add the Scatter trace to the figure
                fig2.add_trace(scatter_trace)

        # plot variant to highlight
        if highlight_var is not None and highlight_var != []:
            highlight_trace = go.Scatter(
                x=subset_var_highlight_df[x_col],
                y=subset_var_highlight_df[y_col],
                mode='markers',
                opacity=1,
                marker=dict(
                    size=7,
                    line=dict(width=3, color=yellow),
                    autocolorscale=True,
                    color=yellow),
                customdata=subset_var_highlight_df[hover_columns],
                hovertemplate="<br>".join(hover_text),
                name="Highlited variants")
            fig2.add_trace(highlight_trace)

    dict_label_axis = {'CADD.phred': 'CADD phred', 'VARITY_R': 'VARITY', 'REVEL': 'REVEL', 'max_spliceAI': 'SpliceAI', 'function_score_final': 'SGE Function Score', 'rna_score': 'RNA score'}
    dict_label_leg = {'clinvar_simple': 'ClinVar', 'consequence': 'Consequence', 'tier_class': 'Function Class', 'Cancer_type_single': 'Cancer Type'}

    # Make it looks cute
    fig2.update_layout(plot_bgcolor=dark_gray,
                       xaxis_title=dict(text=dict_label_axis[x_col], font=dict(color=yel)),
                       yaxis_title=dict(text=dict_label_axis[y_col], font=dict(color=yel)),
                       xaxis=black3dbg,
                       yaxis=black3dbg,
                       paper_bgcolor='rgb(41,41,41)',
                       font_family=font_list[idx_font],
                       font_color=yel,
                       title_font_family=font_list[idx_font],
                       title_text=title,
                       title_y=0.955,
                       title_font_color=yel,
                       title_font_size=18,
                       showlegend=True,
                       height=650,
                       legend=dict(orientation='v', yanchor='bottom', y=0, xanchor='left', x=0, title=dict_label_leg[color_column],
                                   font=dict(color=yel), bgcolor='rgba(41,41,41,0.65)'),
                       modebar=dict(
                           bgcolor=transparent,
                           activecolor=yel,
                           color=yellow)
                       )

    return fig2.update_layout(uirevision=True)


@app.callback(
    Output('dashbio-default-molecule3d', 'modelData'),
    Output('dashbio-default-molecule3d', 'styles'),
    Input('pdb-selector', 'value'),
    Input('vizua_type_3d', 'value'),
    Input(variant_highlight_dropd, 'value'),
)
def update_stucture_based_dropdown(selected_pdb_file, vizu_type, highlight_var):
    data = get_structure_file(selected_pdb_file)
    styles = create_style_3d(
        df, 'average_fs_missense_at_aa_rna', data['atoms'], visualization_type=vizu_type,
        color_element='residue_score', hightlight_vars=highlight_var)
    return data, styles


@app.callback(
    Output(variant_highlight_dropd, 'value'),
    Input('dashbio-default-molecule3d', 'selectedAtomIds'),
    Input('pdb-selector', 'value'),
)
def update_dropdown_based_stucture(atom_ids, selected_pdb_file):
    data = get_structure_file(selected_pdb_file)
    list_var = []  # variants list selected to put in dropdown
    # Get residue index from 60 to 209 to match the structure when VHL
    for elem in data['atoms']:
        if elem['chain'] == 'V':
            elem['residue_index'] = elem['residue_index'] + 60
        else:
            elem['residue_index'] = -1

    if atom_ids is not None and len(atom_ids) > 0:
        last_atom_dict = data['atoms'][atom_ids[-1]]
        subset_df = df.loc[
            (df['protPos'] == last_atom_dict['residue_index']) & (~df['average_fs_missense_at_aa_rna'].isna())]

        if len(list(subset_df['variant_id'])) >= 1:
            list_var = list(subset_df['variant_id'])
    return list_var


@app.callback(
    Output('default-molecule3d-output', 'children'),
    Input('dashbio-default-molecule3d', 'selectedAtomIds'),
    Input('pdb-selector', 'value')
)
def show_selected_residue(atom_ids, selected_pdb_file):
    data = get_structure_file(selected_pdb_file)

    chain_dict = {'H': 'HIF 1A', 'V': 'VHL', 'C': "ELOC", 'B': "ELOB"}
    # Get residue index from 60 to 209 to match the structure when VHL
    for elem in data['atoms']:
        if elem['chain'] == 'V':
            elem['residue_index'] = elem['residue_index'] + 60
        else:
            elem['residue_index'] = -1

    if atom_ids is None or len(atom_ids) == 0:
        return 'Click somewhere on the VHL protein structure to select an amino acid.'

    else:
        last_atom_dict = data['atoms'][atom_ids[-1]]
        # return Only protein / Chain when not VHL
        if (last_atom_dict['chain'] == 'C') or (last_atom_dict['chain'] == 'H') or (last_atom_dict['chain'] == 'B'):
            prot = 'Chain: ', chain_dict[str(last_atom_dict['chain'])],
            phr1 = 'Click somewhere on the VHL protein structure to select an amino acid.'
            return html.Div([html.Br(), html.Div(prot), html.Br(), html.Div(phr1), html.Br()])

    aa_name = 'Reference amino acid: ', last_atom_dict['residue_name'], \
        ', position: ', str(last_atom_dict['residue_index'])
    subset_df = df.loc[(df['protPos'] == last_atom_dict['residue_index']) &
                       (~df['average_fs_missense_at_aa_rna'].isna())]
    return print_var_score_for_selected_residue(subset_df, aa_name)


# Run app
if __name__ == '__main__':
    app.run_server(debug=True)
