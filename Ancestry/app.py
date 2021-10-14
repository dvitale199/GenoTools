####To run this app, you need folloing items####
#coriell.QC.metrics.h5
#style.css
#gp2_2.jpg
#Redlat.png
###Please change the paths
###Please change df_1 as Coriell gwas sum stat
###Change the color of confusion matrix
###There is a known bug with Zoom issue on Manhattan plot, need to be fixed

#upload some libraries
import os
import sys
import subprocess
import numpy as np
import pandas as pd
import pandas.testing as tm
import h5py
import streamlit as st
import matplotlib.pyplot as plt
from matplotlib import style
from matplotlib import cm
import plotly.express as px
import plotly.io as pio
from plotly.subplots import make_subplots
import plotly.graph_objects as go
import seaborn as sns
from PIL import Image
import datetime
import math
from scipy import stats
import statsmodels.formula.api as smf
import statsmodels.api as sm
import statistics
import umap.umap_ as umap
#import umap as umap
from joblib import dump, load
from sklearn.decomposition import PCA
from functools import reduce
import dash_bio as dashbio


st.set_page_config(
    layout = 'wide'
)
# set paths
basedir = '/Users/songy4/Documents'
datadir = f'{basedir}/Ancestry-Estimation-main'

#read corielle.QC.metrics.h5 file
df_qc = pd.read_hdf(f'{datadir}/coriell.QC.metrics.h5', key='QC')
df_ancestry_counts = pd.read_hdf(f'{datadir}/coriell.QC.metrics.h5', key='ancestry_counts')
df_ancestry_labels = pd.read_hdf(f'{datadir}/coriell.QC.metrics.h5', key='ancestry_labels')
df_confusion_matrix = pd.read_hdf(f'{datadir}/coriell.QC.metrics.h5', key='confusion_matrix')
df_new_samples_umap = pd.read_hdf(f'{datadir}/coriell.QC.metrics.h5', key='new_samples_umap')
df_projected_pcs = pd.read_hdf(f'{datadir}/coriell.QC.metrics.h5', key='projected_pcs')
df_ref_pcs = pd.read_hdf(f'{datadir}/coriell.QC.metrics.h5', key='ref_pcs')
df_ref_umap = pd.read_hdf(f'{datadir}/coriell.QC.metrics.h5', key='ref_umap')
df_total_umap = pd.read_hdf(f'{datadir}/coriell.QC.metrics.h5', key='total_umap')

df_1 = pd.read_csv('https://raw.githubusercontent.com/plotly/dash-bio-docs-files/master/manhattan_data.csv')

##Background color
def local_css(file_name):
    with open(file_name) as f:
        st.markdown(f'<style>{f.read()}</style>', unsafe_allow_html=True)
local_css("style.css")
####################### HEAD ##############################################

head_1, head_2, title, head_3, head_4 = st.beta_columns([0.8,1.3, 4, 1, 1])

gp2 = Image.open('gp2_2.jpg')
head_1.image(gp2, width = 120)
redlat = Image.open('Redlat.png')
head_2.image(redlat, width = 180)

with title:
    st.markdown("""
    <style>
    .big-font {
        font-family:IBM Plex Sans; color:#0f557a; font-size:48px !important;
    }
    </style>
    """, unsafe_allow_html=True)
    st.markdown('<p class="big-font">Genetic Report</p>', unsafe_allow_html=True)

with head_3:
    #Study name: get from file name
    a = open('coriell.QC.metrics.h5')
    a_name = a.name
    name, split = a_name.split('.', 1)
#    st.markdown("**Report Full Name**")
#    st.write(a_name)
    st.markdown("**COHORT NAME**")
    st.markdown(name.upper())
    
with head_4:
    def modification_date(filename):
        t = os.path.getmtime(filename)
        return datetime.datetime.fromtimestamp(t)
    date = modification_date(f'{datadir}/coriell.QC.metrics.h5') 
    st.markdown("**REPORT DATE**")  
    st.markdown(date)

########################  SIDE BAR #########################################
st.sidebar.markdown('**Graph selection**', unsafe_allow_html=True)
#prune_selection = df_qc['step'].unique().tolist()
#prune_selection.remove('variant_prune')
selected_metrics = st.sidebar.selectbox(label="Prune selection", options=['All Sample Prune', 'Related Prune', 'Varient Prune'])
#selected_metrics = df_qc.loc[df_qc['step'] == prune_selection]
#selected_metrics = selected_metrics.reset_index()
selected_metrics_1 = st.sidebar.selectbox(label="PCA selection", options=['Reference PCA', 'Projected PCA'])

selected_metrics_2 = st.sidebar.selectbox(label="UMAP selection", options=['Reference UMAP', 'New Sample UMAP', 'Total UMAP'])


    
########################  right column  #########################################
left_column, right_column = st.beta_columns([1.25,1])
with right_column:
    
    #pie_chart for ancestry_counts
    df_ancestry_counts_percent = df_ancestry_counts.copy()
    df_ancestry_counts_percent['percent'] = 100 * df_ancestry_counts_percent['count']  / df_ancestry_counts_percent['count'].sum()
    pie_chart = px.pie(df_ancestry_counts_percent, values=df_ancestry_counts_percent.percent, names = df_ancestry_counts_percent.label)
    pie_chart.update_layout(showlegend=True,
		width=550,
		height=550)
    st.markdown("**Ancestry Distribution**")
    st.write(pie_chart)

        
    #heatmap for confusion matrix
    heatmap = px.imshow(df_confusion_matrix,
                labels=dict(x="Ancestry", y="Ancestry", color="Number"),
                x=df_confusion_matrix.columns.values.tolist(),
                y=df_confusion_matrix.index.tolist()
               )
    heatmap.update_xaxes(side="top")
    st.markdown("**Confusion Matrix**")
    st.write(heatmap)
    
#failed report
    df_6 = df_qc.loc[df_qc['pass'] == False]
    df_6 = df_6.reset_index(drop=True)
    st.markdown("**Failed Prune Steps**")
    st.table(df_6)    

    
########################  left column   #########################################    

with left_column:
#    fig = go.Figure()
    if selected_metrics == 'All Sample Prune':
        #all sample prune
        df_2 = df_qc.query("level == 'sample'")
        df_2['sum'] = df_2.groupby('step')['pruned_count'].transform('sum')
        df_2 = df_2[['step','sum']]
        df_2 = df_2.drop_duplicates(subset=['step', 'sum'], keep='first')

        print(df_2.head())

        #simple bar chart for callrate_prune and sex_prune

        bar_2 = px.bar(df_2, x='step', y='sum', text='sum',
                       hover_data=['step','sum'], 
                       labels={'step':'Pruning Step', 'sum':'Count'}, height=500, width=800)

        #customize figure
        bar_2.update_traces(marker_color='rgb(158,202,225)', marker_line_color='rgb(8,48,107)',
                          marker_line_width=1.5, opacity=0.6)
        st.markdown("**All Sample Prune Count**")
        st.write(bar_2)


    if selected_metrics == 'Related Prune':
        #related prune
        df_3 = df_qc.query("step == 'related_prune'")
        df_3 = df_3[['ancestry', 'pruned_count', 'metric']]

        df_3_related = df_3.query("metric == 'related_count'").reset_index(drop=True)
        df_3_related = df_3_related.rename(columns={'pruned_count': 'related_count'})
        df_3_related = df_3_related.drop('metric', 1)

        df_3_duplicated = df_3.query("metric == 'duplicated_count'").reset_index(drop=True)
        df_3_duplicated = df_3_duplicated.rename(columns={'pruned_count': 'duplicated_count'})
        df_3_duplicated = df_3_duplicated.drop('metric', 1)

        df_4 = pd.merge(df_3_related, df_3_duplicated, on="ancestry")
        df_4.set_index('ancestry', inplace=True)

        bar_3 = go.Figure(data=[
        go.Bar(y=df_4.index, x=df_4['related_count'], orientation='h', name="Related Count", base=0),
        go.Bar(y=df_4.index, x=-df_4['duplicated_count'], orientation='h', name="Duplicated Count", base=0)
        ])

        bar_3.update_layout(
        barmode='stack')

        bar_3.update_yaxes(
            ticktext=df_4.index,
            tickvals=df_4.index
        )
        st.markdown("**Related Prune Count per Ancestry**")
        st.write(bar_3)

    if selected_metrics == 'Varient Prune':
        #vaient prune
        df_5 = df_qc.query("step == 'variant_prune'")
        df_5 = df_5[['ancestry', 'pruned_count', 'metric']]

        df_5_geno = df_5.query("metric == 'geno_removed_count'").reset_index(drop=True)
        df_5_geno = df_5_geno.rename(columns={'pruned_count': 'geno_removed_count'})
        df_5_geno = df_5_geno.drop('metric', 1)

        df_5_mis = df_5.query("metric == 'mis_removed_count'").reset_index(drop=True)
        df_5_mis = df_5_mis.rename(columns={'pruned_count': 'mis_removed_count'})
        df_5_mis = df_5_mis.drop('metric', 1)

        df_5_haplo = df_5.query("metric == 'haplotype_removed_count'").reset_index(drop=True)
        df_5_haplo = df_5_haplo.rename(columns={'pruned_count': 'haplotype_removed_count'})
        df_5_haplo = df_5_haplo.drop('metric', 1)

        df_5_hwe = df_5.query("metric == 'hwe_removed_count'").reset_index(drop=True)
        df_5_hwe = df_5_hwe.rename(columns={'pruned_count': 'hwe_removed_count'})
        df_5_hwe = df_5_hwe.drop('metric', 1)

        df_5_total = df_5.query("metric == 'total_removed_count'").reset_index(drop=True)
        df_5_total = df_5_total.rename(columns={'pruned_count': 'total_removed_count'})
        df_5_total = df_5_total.drop('metric', 1)

        data = [df_5_geno, df_5_mis, df_5_haplo, df_5_hwe, df_5_total]
        df_merged = reduce(lambda left,right: pd.merge(left,right,on=['ancestry'], how='outer'), data)
        df_merged.set_index('ancestry', inplace=True)

        bar_6 = go.Figure(go.Bar(x=df_merged.index, y=df_merged['geno_removed_count'], name='Geno Removed Count'))
        bar_6.add_trace(go.Bar(x=df_merged.index, y=df_merged['mis_removed_count'], name='Mis Removed Count'))
        bar_6.add_trace(go.Bar(x=df_merged.index, y=df_merged['haplotype_removed_count'], name='Haplotype Removed Count'))
        bar_6.add_trace(go.Bar(x=df_merged.index, y=df_merged['hwe_removed_count'], name='Hwe removed Count'))

        bar_6.update_layout(
            xaxis=dict(
                title='Ancestry',
                tickfont_size=14,
            ),
            yaxis=dict(
                title='Count',
                titlefont_size=16,
                tickfont_size=14,
            )
        )

        bar_6.update_layout(barmode='stack', xaxis={'categoryorder':'total descending'},width=800,height=600)
        st.markdown("**Varient Prune Count per Ancestry**")
        st.write(bar_6)
        
#PCA
    if selected_metrics_1 == 'Reference PCA':
        X = df_ref_pcs[['PC1', 'PC2']]
        pca = PCA(n_components=2)
        components = pca.fit_transform(X)
        pcs_1 = px.scatter(components, x=0, y=1, color=df_ref_pcs['label'],
                          labels={'0': 'Principal Component 1', '1': 'Principal Component 2'})
        st.markdown("**Reference PCA**")
        st.write(pcs_1)
        
    if selected_metrics_1 == 'Projected PCA':
        X = df_projected_pcs[['PC1', 'PC2']]
        pca = PCA(n_components=2)
        components = pca.fit_transform(X)
        pcs_2 = px.scatter(components, x=0, y=1, color=df_projected_pcs['label'],
                          labels={'0': 'Principal Component 1', '1': 'Principal Component 2'})
        st.markdown("**Projected PCA**")
        st.write(pcs_2)
            
        
#UMAP
    if selected_metrics_2 == 'Reference UMAP':
        N = 3
        features = df_ref_umap.iloc[:, :N]
        #umap_3d = umap.UMAP()
        umap_3d = umap.UMAP(n_components=3, init='random', random_state=0)
        proj_3d = umap_3d.fit_transform(features)
        umap_1 = px.scatter_3d(
            proj_3d, x=0, y=1, z=2,
            color=df_ref_umap.label, 
            labels={'0': 'UMAP_1', '1': 'UMAP_2', '2': 'UMAP_3'}
        )
        umap_1.update_layout(showlegend=True,
		width=550,
		height=550)
        st.markdown("**Reference UMAP**")
        st.write(umap_1)
        
    if selected_metrics_2 == 'New Sample UMAP':
        N = 3
        features = df_new_samples_umap.iloc[:, :N]
        #umap_3d = umap.UMAP()
        umap_3d = umap.UMAP(n_components=3, init='random', random_state=0)
        proj_3d = umap_3d.fit_transform(features)
        umap_2 = px.scatter_3d(
            proj_3d, x=0, y=1, z=2,
            color=df_new_samples_umap.label, 
            labels={'0': 'UMAP_1', '1': 'UMAP_2', '2': 'UMAP_3'}
        )
        umap_2.update_layout(showlegend=True,
		width=550,
		height=550)
        st.markdown("**New Sample UMAP**")
        st.write(umap_2)
        
    if selected_metrics_2 == 'Total UMAP':
        N = 3
        features = df_total_umap.iloc[:, :N]
        #umap_3d = umap.UMAP()
        umap_3d = umap.UMAP(n_components=3, init='random', random_state=0)
        proj_3d = umap_3d.fit_transform(features)
        umap_3 = px.scatter_3d(
            proj_3d, x=0, y=1, z=2,
            color=df_total_umap.label, 
            labels={'0': 'UMAP_1', '1': 'UMAP_2', '2': 'UMAP_3'}
        )
        umap_3.update_layout(showlegend=True,
		width=550,
		height=550)
        st.markdown("**Total UMAP**")
        st.write(umap_3)

        
########################  SIDE BAR #########################################

#df_1['CHR_NUM'] = 'Chr'+ df_1['CHR'].astype(str)
#chromosome_selected = st.multiselect('CHR', chromosome)
#if len(chromosome_selected) > 0:
#    df_cf = df_1.loc[df_1['CHR'].isin(countries_selected)]
#else:
#    df_cf = df


threshold = st.slider('Threshold value', min_value = 1, max_value = 10, value =1)

             

#manhattan plot
    
DATASET = df_1.groupby('CHR').apply(lambda u: u.head(50))
DATASET = DATASET.droplevel('CHR').reset_index(drop=True)

figure = dashbio.ManhattanPlot(
    dataframe = df_1,
    highlight_color = "#ff3860",
    genomewideline_value = threshold,
    genomewideline_width = 2,
    genomewideline_color = "#00d1b2",
    suggestiveline_value = False,
    showgrid = False,
    title = None,
    xlabel = "chromosome"
)
figure.update_layout(showlegend=True,
    width=1200,
    height=500)
st.markdown("**Manhattan Plot**")
st.write(figure)