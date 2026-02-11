# Methodological considerations in the analysis of guard cell starch metabolism

Statistics and figure generation for the correspondence: *Methodological considerations in the analysis of guard cell starch metabolism*.

## Structure

- `code/` - R and Python scripts for plot generation and formatting
- `assets/plots/` - Generated figures (PDF and PNG)
- `assets/legends/` - Figure legend text (markdown)
- `data/` - Raw data files

## Figures

- **Figure 1**: GC starch dynamics and stomatal aperture in response to blue light
- **Figure 2**: GC viability in different stomatal opening buffers and leaf starch amounts
- **Figure S1**: Stomata width analysis

Each figure is output as: original PDF, A4 PDF, A4 PDF with legend, and 600 DPI PNG.

## Usage

```r
source("code/combined_plots.r")           # Generate plots from R
```
```bash
python code/create_a4_versions.py         # A4 PDFs without legends
python code/create_a4_with_legends.py     # A4 PDFs with legends
python code/create_high_dpi_png.py        # 600 DPI PNGs
python code/crop_annotated_objects.py     # Crop annotated objects from images
```

## Requirements

**R**: ggplot2, dplyr, tidyr, gridExtra, grid

**Python**: `pip install pikepdf reportlab pdf2image pillow`
