"""
Convert PDF plots to high-DPI PNG files without legends.
Uses pdf2image to convert original plot PDFs to PNG at 300 DPI.
"""

from pathlib import Path
from pdf2image import convert_from_path
import os


def pdf_to_png(pdf_path, output_path, dpi=300):
    """
    Convert a PDF to high-DPI PNG.

    Args:
        pdf_path: Path to the PDF file
        output_path: Path to save the PNG file
        dpi: Resolution in dots per inch (default: 300)
    """
    # Convert PDF to images (returns a list of PIL Image objects)
    images = convert_from_path(pdf_path, dpi=dpi)

    # Save the first page as PNG
    if images:
        images[0].save(output_path, 'PNG')
        print(f"Created high-DPI PNG: {output_path}")
        print(f"  Resolution: {dpi} DPI")
        print(f"  Size: {images[0].width} x {images[0].height} pixels")
    else:
        print(f"Error: No images found in {pdf_path}")


def main():
    plots_dir = Path('assets/plots')
    dpi = 600  # Very high-quality DPI for publication

    print("=" * 80)
    print(f"CREATING HIGH-DPI PNG FILES (DPI: {dpi})")
    print("=" * 80)
    print()

    # Figure 1
    pdf_to_png(
        plots_dir / 'figure_1.pdf',
        plots_dir / 'figure_1.png',
        dpi=dpi
    )
    print()

    # Figure 2
    pdf_to_png(
        plots_dir / 'figure_2.pdf',
        plots_dir / 'figure_2.png',
        dpi=dpi
    )
    print()

    # Supplementary Figure S1
    pdf_to_png(
        plots_dir / 'figure_S1_stomata_width.pdf',
        plots_dir / 'figure_S1_stomata_width.png',
        dpi=dpi
    )
    print()

    print("=" * 80)
    print("Done! Created high-DPI PNG files")
    print("=" * 80)


if __name__ == '__main__':
    main()
