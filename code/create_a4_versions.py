"""
Create A4-sized versions of the original combined plot PDFs.
Embeds the original vectorized PDFs into A4 pages with proper scaling and margins.
Uses pikepdf for reliable PDF manipulation.
"""

from pathlib import Path
import pikepdf
from pikepdf import Pdf, Page, Rectangle

def embed_pdf_in_a4(input_pdf_path, output_pdf_path, margin_inches=0.5):
    """
    Embed a PDF into an A4-sized page with margins, maintaining aspect ratio.

    Args:
        input_pdf_path: Path to the original PDF
        output_pdf_path: Path to save the A4 version
        margin_inches: Margin size in inches (default: 0.5)
    """
    # A4 dimensions in points (1 inch = 72 points)
    a4_width = 8.27 * 72  # 595.27 points
    a4_height = 11.69 * 72  # 841.89 points
    margin_pts = margin_inches * 72

    # Open the original PDF
    src = Pdf.open(input_pdf_path)
    dst = Pdf.new()

    # Get first page to calculate dimensions
    first_page = src.pages[0]
    mediabox = first_page.mediabox
    orig_width = float(mediabox[2] - mediabox[0])
    orig_height = float(mediabox[3] - mediabox[1])

    # Calculate available space after margins
    available_width = a4_width - 2 * margin_pts
    available_height = a4_height - 2 * margin_pts

    # Calculate scaling factor to fit while maintaining aspect ratio
    scale_factor = min(available_width / orig_width, available_height / orig_height)

    # Calculate final dimensions
    final_width = orig_width * scale_factor
    final_height = orig_height * scale_factor

    # Calculate position to center the content
    x_offset = (a4_width - final_width) / 2
    y_offset = (a4_height - final_height) / 2

    # Process each page
    for page in src.pages:
        # Set the page's MediaBox to A4 size FIRST
        page.mediabox = Rectangle(0, 0, a4_width, a4_height)

        # Get the content stream and prepend transformation matrix
        # PDF transformation matrix format: a b c d e f cm
        # For scaling and translation: scale 0 0 scale tx ty cm
        transform_cmd = f"q {scale_factor} 0 0 {scale_factor} {x_offset} {y_offset} cm\n"
        restore_cmd = " Q"

        # Get existing content
        contents = page.obj.Contents
        if isinstance(contents, pikepdf.Array):
            # Multiple content streams - get first one
            old_stream = contents[0].read_bytes()
        else:
            # Single content stream
            old_stream = contents.read_bytes()

        # Create new content with transformation
        new_content = transform_cmd.encode('latin-1') + old_stream + restore_cmd.encode('latin-1')

        # Replace the content stream
        page.obj.Contents = src.make_stream(new_content)

        # Add page to destination
        dst.pages.append(page)

    # Save the output PDF
    dst.save(output_pdf_path)

    print(f"Created A4 version: {output_pdf_path}")
    print(f"  Original: {orig_width/72:.2f} x {orig_height/72:.2f} inches")
    print(f"  Scale factor: {scale_factor:.3f}")
    print(f"  Final size: {final_width/72:.2f} x {final_height/72:.2f} inches (on A4 page)")


def main():
    # Define paths
    plots_dir = Path('assets/plots')

    # List of figures to process
    figures = [
        'figure_1.pdf',
        'figure_2.pdf',
        'figure_S1_stomata_width.pdf'
    ]

    print("=" * 80)
    print("CREATING A4-SIZED VERSIONS (VECTORIZED)")
    print("=" * 80)
    print()

    for figure in figures:
        input_path = plots_dir / figure
        output_path = plots_dir / figure.replace('.pdf', '_A4.pdf')

        if input_path.exists():
            embed_pdf_in_a4(input_path, output_path, margin_inches=0.5)
            print()
        else:
            print(f"Warning: {input_path} not found, skipping...")

    print("=" * 80)
    print("Done! All A4 versions created with _A4.pdf suffix")
    print("All files remain fully vectorized (no rasterization)")
    print("=" * 80)


if __name__ == '__main__':
    main()
