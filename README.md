# Daloso Project

Scientific figure generation pipeline for guard cell starch dynamics and stomatal aperture research.

## Project Structure

```
Daloso/
├── code/                           # Analysis and figure generation scripts
│   ├── combined_plots.r           # R script to generate original plots
│   ├── create_a4_versions.py      # Create A4 PDFs without legends
│   ├── create_a4_with_legends.py  # Create A4 PDFs with Arial legends
│   └── create_high_dpi_png.py     # Generate 600 DPI PNG files
├── assets/
│   ├── legends/                   # Figure legend text (markdown format)
│   │   ├── figure_1_legend.md
│   │   └── figure_2_legend.md
│   └── plots/                     # Generated figures (PDFs and PNGs)
└── data/                          # Raw data files
```

## Generated Figures

### Main Figures
- **Figure 1**: GC starch dynamics and stomatal aperture in response to blue light
- **Figure 2**: GC viability in different stomatal opening buffers and leaf starch amounts

### Supplementary Figures
- **Figure S1**: Stomata width analysis

## Output Formats

For each figure, multiple formats are generated:

1. **Original PDF** (`figure_X.pdf`) - 14 × 18 inches vectorized plot
2. **A4 PDF** (`figure_X_A4.pdf`) - A4-sized without legend
3. **A4 PDF with legend** (`figure_X_A4_legend.pdf`) - A4-sized with formatted legend (9pt Arial)
4. **High-DPI PNG** (`figure_X.png`) - 600 DPI raster image (8400 × 10800 pixels)

## Usage

### Generate all figures from R

```r
source("code/combined_plots.r")
```

### Create A4 versions without legends

```bash
python code/create_a4_versions.py
```

### Create A4 versions with legends

```bash
python code/create_a4_with_legends.py
```

### Generate high-resolution PNG files

```bash
python code/create_high_dpi_png.py
```

## Requirements

### R Dependencies
- ggplot2
- dplyr
- tidyr
- gridExtra
- grid

### Python Dependencies
- pikepdf
- reportlab
- pdf2image

Install Python dependencies:
```bash
pip install pikepdf reportlab pdf2image
```

## Legend Formatting

Legend text files use markdown formatting:
- `**bold text**` → Bold
- `*italic text*` → Italic
- `\*` → Literal asterisk (for statistical notation like `*, p < 0.05`)

Legends are automatically converted to HTML and rendered in Arial font with justified alignment.

## Technical Details

### A4 PDF Generation
- A4 dimensions: 8.27 × 11.69 inches
- Plot margins: 0.1 inches (minimized for maximum plot size)
- Legend margins: 0.25 inches
- Legend font: Arial 9pt (with bold/italic support)
- Spacing between plot and legend: 0.15 inches

### High-DPI PNG Settings
- Resolution: 600 DPI
- Pixel dimensions: 8400 × 10800 pixels
- Format: PNG (lossless compression)
- Source: Original vectorized PDFs

## Notes

- All PDFs maintain vector graphics (no rasterization)
- Plot fonts are preserved from R-generated originals
- Legend fonts use Arial with proper bold/italic rendering
- Font resources are properly merged to support formatted text
