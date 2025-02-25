library(DiagrammeR)

grViz("
  digraph G {
    rankdir=TB
    graph [bgcolor = '#FDFDFD']
    edge [style=invis]
    node [fontname = 'helvetica', width = 3, height = 1, fontsize=12, fixedsize=true]

    subgraph cluster_data {
      label = 'Data Selection & Preprocessing'
      style=dashed
      color= '#625a5a'
      fontname = 'helvetica-bold'

      node [shape=rectangle, style=filled, color=black, fillcolor = '#91cf60', fontcolor = '#2a6496']
      tcga_samples [label = 'TCGA colon adenocarcinoma\nsamples (n = 524)']
      filter_primary [label = 'Filter out normal and metastasis\nsamples (n = 481)']
      filter_clinical [label = 'Filter out samples\nwith missing clinical data\n(n = 475: 225F, 250M)']
      normalization [label = 'Normalization\nand filtering']

      tcga_samples -> filter_primary -> filter_clinical -> normalization
    }

    subgraph cluster_dea {
      label = 'Differential Expression Analysis'
      style=dashed
      color= '#625a5a'
      fontname = 'helvetica-bold'

      node [shape=rectangle, style=filled, color=black, fillcolor = '#fee08b', fontcolor = '#2a6496']
      dea_analysis [label = '327 DEGs\n(122 Down, 205 Up)']

      normalization -> dea_analysis
    }

    # Functional Enrichment Analysis and Survival Analysis directly under DEA
    subgraph cluster_analysis {
      label = ''
      style=dashed
      color= '#625a5a'
      fontname = 'helvetica-bold'

      node [shape=rectangle, style=filled, color=black, fillcolor = '#d9f5d3', fontcolor = '#2a6496']
      enrichment [label = 'Functional enrichment\nanalysis']

      node [shape=rectangle, style=filled, color=black, fillcolor = '#f9d3d3', fontcolor = '#2a6496']
      survival [label = 'Survival analysis\n(CoxPH & KM)']

      dea_analysis -> enrichment
      dea_analysis -> survival

      { rank = same; enrichment; survival; }  # Aligning enrichment and survival horizontally
    }

    subgraph cluster_prognostic {
      label = ''
      style=dashed
      color= '#625a5a'
      fontname = 'helvetica-bold'

      node [shape=rectangle, style=filled, color=black, fillcolor = '#ffcccb', fontcolor = '#2a6496']
      prognostic_markers [label = 'Identification of clinically significant\nprognostic Markers']

      survival -> prognostic_markers
    }

    # Aesthetic edge adjustments
    node [shape=rectangle, fontsize=11]
    edge [color='#808080', fontname='Helvetica', style=solid, arrowhead='vee']

    # Connecting nodes with arrows
    tcga_samples -> filter_primary [color='#2a6496', arrowsize=1.2]
    filter_primary -> filter_clinical [color='#2a6496', arrowsize=1.2]
    filter_clinical -> normalization [color='#2a6496', arrowsize=1.2]
    normalization -> dea_analysis [color='#ffb84d', arrowsize=1.2]
    dea_analysis -> enrichment [color='#72c47e', arrowsize=1.2]
    dea_analysis -> survival [color='#72c47e', arrowsize=1.2]
    survival -> prognostic_markers [color='#ff6f69', arrowsize=1.2]
  }
")

