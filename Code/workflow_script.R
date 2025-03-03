library(DiagrammeR)

grViz("
  digraph G {
    rankdir=TB
    graph [bgcolor = '#FDFDFD']
    edge [style = invis]
    node [fontname = 'helvetica', width = 3, height = 1, fontsize=12, fixedsize=true]

    subgraph cluster_data {
      label = 'Data Selection & Preprocessing'
      style=dashed
      color= '#625a5a'
      fontname = 'helvetica-bold'

      node [shape=rectangle, style=filled, color=black, fillcolor = '#91cf60', fontcolor = '#2a6496']
      tcga_samples [label = 'TCGA-COAD RNA-seq\n(n = 524)']
      filter_primary [label = 'Primary tumors\n(n = 474)']
      filter_unnecessary [label = 'Filtering unnecessary samples\n(n = 465: 221F, 244M)']
      normalization [label = 'Normalization\nand filtering']

      tcga_samples -> filter_primary -> filter_unnecessary -> normalization
    }

    subgraph cluster_dea {
      label = 'Differential Expression Analysis'
      style=dashed
      color= '#625a5a'
      fontname = 'helvetica-bold'

      node [shape=rectangle, style=filled, color=black, fillcolor = '#fee08b', fontcolor = '#2a6496']
      dea_analysis [label = '325 DEGs (128 up, 197 down)\n(FDR < 0.01, |logFC| > 1)']

      normalization -> dea_analysis
    }

    # Functional Enrichment Analysis Cluster
    subgraph cluster_enrichment {
      label = 'Functional Enrichment Analysis'
      style=dashed
      color= '#625a5a'
      fontname = 'helvetica-bold'

      node [shape=rectangle, style=filled, color=black, fillcolor = '#d9f5d3', fontcolor = '#2a6496']
      enrichment [label = 'KEGG Pathways & GO']

      dea_analysis -> enrichment
    }

    # Survival Analysis Cluster
    subgraph cluster_survival {
      label = 'Survival Analysis'
      style=dashed
      color= '#625a5a'
      fontname = 'helvetica-bold'

      node [shape=rectangle, style=filled, color=black, fillcolor = '#ffcccb', fontcolor = '#2a6496']
      survival [label = 'Cox Proportional Hazards Analysis']

      dea_analysis -> survival
    }

    # Male Cohort Analysis
    subgraph cluster_male {
      label = 'Male Cohort'
      style=dashed
      color= '#625a5a'
      fontname = 'helvetica-bold'

      node [shape=rectangle, style=filled, color=black, fillcolor = '#ff9999', fontcolor = '#2a6496']
      uva_male [label = 'UVA Cox Analysis']
      mva_male [label = 'MVA Cox Analysis\nwith\nBackward Selection']
      km_male [label = 'Kaplan-Meier Analysis']

      survival -> uva_male -> mva_male -> km_male
    }

    # Female Cohort Analysis
    subgraph cluster_female {
      label = 'Female Cohort'
      style=dashed
      color= '#625a5a'
      fontname = 'helvetica-bold'

      node [shape=rectangle, style=filled, color=black, fillcolor = '#ffb3e6', fontcolor = '#2a6496']
      uva_female [label = 'UVA Cox Analysis']
      mva_female [label = 'MVA Cox Analysis\nwith\nBackward Selection']
      km_female [label = 'Kaplan-Meier Analysis']

      survival -> uva_female -> mva_female -> km_female
    }

    # Aesthetic edge adjustments
    node [shape=rectangle, fontsize=11]
    edge [color='#808080', fontname='Helvetica', style=solid, arrowhead='vee']

    # Connecting nodes with arrows
    tcga_samples -> filter_primary [color='#2a6496', arrowsize=1.2]
    filter_primary -> filter_unnecessary [color='#2a6496', arrowsize=1.2]
    filter_unnecessary -> normalization [color='#2a6496', arrowsize=1.2]
    normalization -> dea_analysis [color='#ffb84d', arrowsize=1.2]
    dea_analysis -> enrichment [color='#72c47e', arrowsize=1.2]
    dea_analysis -> survival [color='#72c47e', arrowsize=1.2]
    survival -> uva_male [color='#ff6f69', arrowsize=1.2]
    survival -> uva_female [color='#ff6f69', arrowsize=1.2]
    uva_male -> mva_male [color='#ff6f69', arrowsize=1.2]
    mva_male -> km_male [color='#ff6f69', arrowsize=1.2]
    uva_female -> mva_female [color='#ff6f69', arrowsize=1.2]
    mva_female -> km_female [color='#ff6f69', arrowsize=1.2]
  }
")
