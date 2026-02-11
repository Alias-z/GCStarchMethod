"""
Create A4-sized versions of plots with legends below.
Embeds the original vectorized PDFs and adds formatted legend text.
"""

from pathlib import Path
import pikepdf
from pikepdf import Pdf, Rectangle
from reportlab.pdfgen import canvas
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.platypus import Paragraph
from reportlab.lib.units import inch
from reportlab.lib.enums import TA_JUSTIFY
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
import io
import os


# Register Arial font if available, otherwise fall back to Helvetica
def register_arial():
    """Register Arial font from Windows system fonts."""
    try:
        # Try to find Arial font on Windows
        windows_fonts = r'C:\Windows\Fonts'
        arial_path = os.path.join(windows_fonts, 'arial.ttf')
        arial_bold_path = os.path.join(windows_fonts, 'arialbd.ttf')
        arial_italic_path = os.path.join(windows_fonts, 'ariali.ttf')
        arial_bold_italic_path = os.path.join(windows_fonts, 'arialbi.ttf')

        if os.path.exists(arial_path):
            pdfmetrics.registerFont(TTFont('Arial', arial_path))
            if os.path.exists(arial_bold_path):
                pdfmetrics.registerFont(TTFont('Arial-Bold', arial_bold_path))
            if os.path.exists(arial_italic_path):
                pdfmetrics.registerFont(TTFont('Arial-Italic', arial_italic_path))
            if os.path.exists(arial_bold_italic_path):
                pdfmetrics.registerFont(TTFont('Arial-BoldItalic', arial_bold_italic_path))

            # CRITICAL: Register the font family so <b> and <i> tags work
            from reportlab.pdfbase.pdfmetrics import registerFontFamily
            registerFontFamily('Arial',
                             normal='Arial',
                             bold='Arial-Bold',
                             italic='Arial-Italic',
                             boldItalic='Arial-BoldItalic')
            return 'Arial'
        else:
            return 'Helvetica'  # Fallback to Helvetica (standard PDF font, very similar to Arial)
    except Exception as e:
        print(f"Warning: Could not register Arial font: {e}")
        return 'Helvetica'


# Register the font at module load
LEGEND_FONT = register_arial()


def read_legend(legend_path):
    """Read and clean legend text from markdown file."""
    with open(legend_path, 'r', encoding='utf-8') as f:
        text = f.read().strip()

    import re

    # FIRST: Handle escaped asterisks (convert \* to a placeholder)
    # This preserves literal asterisks in statistical notation
    ASTERISK_PLACEHOLDER = '<!ASTERISK!>'
    text = text.replace('\\*', ASTERISK_PLACEHOLDER)

    # THEN: Convert markdown formatting to HTML
    # Handle bold+italic patterns (***text***)
    text = re.sub(r'\*\*\*(.*?)\*\*\*', r'<b><i>\1</i></b>', text)

    # Handle bold (**text**)
    text = re.sub(r'\*\*(.*?)\*\*', r'<b>\1</b>', text)

    # Handle italic (single asterisk *text*)
    text = re.sub(r'\*([^*]+?)\*', r'<i>\1</i>', text)

    # FINALLY: Restore literal asterisks from placeholders
    text = text.replace(ASTERISK_PLACEHOLDER, '*')

    # Escape ampersands first (but not in HTML entities)
    text = re.sub(r'&(?!amp;|lt;|gt;)', '&amp;', text)

    # Escape < and > except in our HTML tags
    def escape_except_tags(text):
        result = []
        i = 0
        while i < len(text):
            if text[i:i+3] in ['<b>', '<i>']:
                result.append(text[i:i+3])
                i += 3
            elif text[i:i+4] in ['</b>', '</i>']:
                result.append(text[i:i+4])
                i += 4
            elif text[i] == '<':
                result.append('&lt;')
                i += 1
            elif text[i] == '>':
                result.append('&gt;')
                i += 1
            else:
                result.append(text[i])
                i += 1
        return ''.join(result)

    text = escape_except_tags(text)

    return text


def create_legend_pdf(legend_text, width, height, font_size=10, margin_pts=None):
    """Create a PDF with formatted legend text."""
    if margin_pts is None:
        margin_pts = 0.4 * 72  # Default 0.4 inch margin

    packet = io.BytesIO()
    c = canvas.Canvas(packet, pagesize=(width, height))

    # Create custom paragraph style for justified text
    styles = getSampleStyleSheet()
    legend_style = ParagraphStyle(
        'Legend',
        parent=styles['Normal'],
        fontSize=font_size,
        leading=font_size * 1.2,  # Line spacing
        alignment=TA_JUSTIFY,
        fontName=LEGEND_FONT,  # Use Arial if available, Helvetica otherwise
        textColor='black',
        embeddedHyphenation=0
    )

    # Create paragraph
    para = Paragraph(legend_text, legend_style)

    # Calculate available width (with margins)
    available_width = width - 2 * margin_pts

    # Get required height for this text
    w, h = para.wrap(available_width, height)

    # Draw the paragraph at the top of the space
    para.drawOn(c, margin_pts, height - h - margin_pts)

    c.save()
    packet.seek(0)
    return packet, h + 2 * margin_pts  # Return actual height used


def embed_pdf_with_legend(input_pdf_path, legend_path, output_pdf_path,
                          font_size=10, margin_inches=0.25, legend_margin_inches=0.4):
    """
    Embed a PDF and legend into an A4-sized page.

    Args:
        input_pdf_path: Path to the original plot PDF
        legend_path: Path to the legend markdown file
        output_pdf_path: Path to save the A4 version
        font_size: Font size for legend text (default: 10)
        margin_inches: Margin size for plot in inches (default: 0.25)
        legend_margin_inches: Margin size for legend text in inches (default: 0.4)
    """
    # A4 dimensions in points (1 inch = 72 points)
    a4_width = 8.27 * 72  # 595.27 points
    a4_height = 11.69 * 72  # 841.89 points
    margin_pts = margin_inches * 72
    legend_margin_pts = legend_margin_inches * 72

    # Read legend
    legend_text = read_legend(legend_path)

    # Create legend PDF to measure height needed (using legend margins for text)
    legend_pdf_buffer, legend_height = create_legend_pdf(
        legend_text, a4_width, a4_height, font_size, legend_margin_pts
    )

    # Open the original plot PDF
    src = Pdf.open(input_pdf_path)

    # Get first page to calculate dimensions
    first_page = src.pages[0]
    mediabox = first_page.mediabox
    orig_width = float(mediabox[2] - mediabox[0])
    orig_height = float(mediabox[3] - mediabox[1])

    # Calculate available space for plot (leaving room for legend at bottom)
    available_width = a4_width - 2 * margin_pts
    spacing = 0.15 * 72  # Small spacing between plot and legend
    available_height = a4_height - legend_height - 2 * margin_pts - spacing

    # Calculate scaling factor to fit while maintaining aspect ratio
    scale_factor = min(available_width / orig_width, available_height / orig_height)

    # Calculate final plot dimensions
    final_width = orig_width * scale_factor
    final_height = orig_height * scale_factor

    # Calculate position to center the plot horizontally, align to top
    x_offset = (a4_width - final_width) / 2
    y_offset = a4_height - final_height - margin_pts  # Align to top with margin

    # Create output PDF
    dst = Pdf.new()

    # Process the plot page
    for page in src.pages:
        # Set the page's MediaBox to A4 size
        page.mediabox = Rectangle(0, 0, a4_width, a4_height)

        # Get the content stream and prepend transformation matrix
        transform_cmd = f"q {scale_factor} 0 0 {scale_factor} {x_offset} {y_offset} cm\n"
        restore_cmd = " Q"

        # Get existing content
        contents = page.obj.Contents
        if isinstance(contents, pikepdf.Array):
            old_stream = contents[0].read_bytes()
        else:
            old_stream = contents.read_bytes()

        # Create new content with transformation
        new_content = transform_cmd.encode('latin-1') + old_stream + restore_cmd.encode('latin-1')

        # Replace the content stream
        page.obj.Contents = src.make_stream(new_content)

        # Add page to destination
        dst.pages.append(page)

    # Now add the legend at the bottom
    # Create legend PDF again for actual use
    legend_pdf_buffer, _ = create_legend_pdf(legend_text, a4_width, legend_height, font_size, legend_margin_pts)
    legend_pdf = Pdf.open(legend_pdf_buffer)

    # Merge legend onto the page
    final_page = dst.pages[0]
    legend_page = legend_pdf.pages[0]

    # CRITICAL: Merge the font resources from legend page to final page
    # This ensures bold/italic fonts are available WITHOUT overwriting plot fonts
    if '/Resources' in legend_page.obj:
        legend_resources = legend_page.obj.Resources
        if '/Resources' not in final_page.obj:
            final_page.obj.Resources = pikepdf.Dictionary()

        final_resources = final_page.obj.Resources

        # Merge Font resources - ADD to existing fonts, don't replace
        if '/Font' in legend_resources:
            if '/Font' not in final_resources:
                final_resources.Font = pikepdf.Dictionary()

            # Copy each font from legend to final page using copy_foreign
            # Only add fonts that don't already exist (preserve plot fonts)
            for font_name, font_obj in legend_resources.Font.items():
                if font_name not in final_resources.Font:
                    final_resources.Font[font_name] = dst.copy_foreign(font_obj)

    # Position legend at bottom
    legend_y_offset = margin_pts

    # Add transformation to position legend
    legend_transform = f"q 1 0 0 1 0 {legend_y_offset} cm\n"
    legend_restore = " Q"

    # Get legend content
    legend_contents = legend_page.obj.Contents
    if isinstance(legend_contents, pikepdf.Array):
        legend_stream = legend_contents[0].read_bytes()
    else:
        legend_stream = legend_contents.read_bytes()

    # Combine plot and legend content
    final_contents = final_page.obj.Contents
    if isinstance(final_contents, pikepdf.Array):
        plot_stream = final_contents[0].read_bytes()
    else:
        plot_stream = final_contents.read_bytes()

    combined_content = (plot_stream + b"\n" +
                       legend_transform.encode('latin-1') +
                       legend_stream +
                       legend_restore.encode('latin-1'))

    final_page.obj.Contents = dst.make_stream(combined_content)

    # Save the output PDF
    dst.save(output_pdf_path)

    print(f"Created A4 version with legend: {output_pdf_path}")
    print(f"  Original plot: {orig_width/72:.2f} x {orig_height/72:.2f} inches")
    print(f"  Scale factor: {scale_factor:.3f}")
    print(f"  Final plot size: {final_width/72:.2f} x {final_height/72:.2f} inches")
    print(f"  Legend height: {legend_height/72:.2f} inches")
    print(f"  Font size: {font_size}pt")


def main():
    plots_dir = Path('assets/plots')
    legends_dir = Path('assets/legends')

    # Test both font sizes
    for font_size in [9, 12]:
        print("=" * 80)
        print(f"CREATING A4 VERSIONS WITH LEGENDS (Font size: {font_size}pt)")
        print("=" * 80)
        print()

        # Figure 1
        embed_pdf_with_legend(
            plots_dir / 'figure_1.pdf',
            legends_dir / 'figure_1_legend.md',
            plots_dir / f'figure_1_A4_legend_{font_size}pt.pdf',
            font_size=font_size,
            margin_inches=0.25,  # Smaller margin for plot
            legend_margin_inches=0.4  # Reasonable margin for legend text
        )
        print()

        # Figure 2
        embed_pdf_with_legend(
            plots_dir / 'figure_2.pdf',
            legends_dir / 'figure_2_legend.md',
            plots_dir / f'figure_2_A4_legend_{font_size}pt.pdf',
            font_size=font_size,
            margin_inches=0.25,
            legend_margin_inches=0.4
        )
        print()

    print("=" * 80)
    print("Done! Created versions with 9pt and 12pt fonts")
    print("Compare both to choose the best option")
    print("=" * 80)


if __name__ == '__main__':
    main()
