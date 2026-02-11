"""
Script to crop annotated objects from ISAT annotation files.
Each object is first made square, then 25% padding is added before cropping.
"""

import json
import os
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont
import numpy as np


def add_scale_bar(img, scale_bar_um=10, pixels_per_um=8):
    """
    Add a scale bar to the bottom right corner of the image.

    Args:
        img: PIL Image object
        scale_bar_um: Length of scale bar in micrometers (default: 10)
        pixels_per_um: Pixels per micrometer (default: 8)

    Returns:
        PIL Image with scale bar added
    """
    # Create a copy to avoid modifying original
    img_with_bar = img.copy()
    draw = ImageDraw.Draw(img_with_bar)

    # Calculate scale bar dimensions
    bar_length_pixels = int(scale_bar_um * pixels_per_um)
    bar_height = 4  # Height of the scale bar in pixels
    margin = 10  # Margin from edges

    # Position: bottom right
    img_width, img_height = img.size
    x_end = img_width - margin
    x_start = x_end - bar_length_pixels
    y_bar = img_height - margin - bar_height - 15  # Leave space for text

    # Draw white background rectangle for better visibility
    bg_padding = 5
    draw.rectangle(
        [x_start - bg_padding, y_bar - bg_padding,
         x_end + bg_padding, img_height - margin + bg_padding],
        fill='white',
        outline='black',
        width=1
    )

    # Draw the scale bar
    draw.rectangle(
        [x_start, y_bar, x_end, y_bar + bar_height],
        fill='black'
    )

    # Add text label
    text = f"{scale_bar_um} μm"

    # Try to use a font, fall back to default if not available
    try:
        font = ImageFont.truetype("arial.ttf", 12)
    except:
        font = ImageFont.load_default()

    # Get text size and center it above the bar
    try:
        bbox = draw.textbbox((0, 0), text, font=font)
        text_width = bbox[2] - bbox[0]
        text_height = bbox[3] - bbox[1]
    except:
        text_width = len(text) * 6
        text_height = 10

    text_x = x_start + (bar_length_pixels - text_width) // 2
    text_y = y_bar + bar_height + 3

    draw.text((text_x, text_y), text, fill='black', font=font)

    return img_with_bar


def crop_objects_from_annotation(json_path, output_dir, padding_percent=25, target_category="stoma"):
    """
    Crop objects from an image based on ISAT annotation file.
    First makes the bounding box square, then applies padding.

    Args:
        json_path: Path to the JSON annotation file
        output_dir: Directory to save cropped images
        padding_percent: Percentage of padding to add around square bounding box (default: 25)
        target_category: Only crop objects with this category (default: "stoma")
    """
    # Read annotation file
    with open(json_path, 'r') as f:
        data = json.load(f)

    # Get image info
    info = data['info']
    image_name = info['name']

    # Construct the image path (same directory as JSON, but with PNG extension)
    json_dir = os.path.dirname(json_path)
    base_name = os.path.splitext(os.path.basename(json_path))[0]
    image_path = os.path.join(json_dir, f"{base_name}.png")

    # Check if image exists
    if not os.path.exists(image_path):
        print(f"Warning: Image not found at {image_path}")
        return

    # Load the image
    img = Image.open(image_path)
    img_width, img_height = img.size

    # Process each object (filter by target category)
    objects = data['objects']
    filtered_objects = [obj for obj in objects if obj['category'] == target_category]

    if not filtered_objects:
        print(f"No '{target_category}' objects found in {image_path}")
        return

    print(f"Processing {image_path} with {len(filtered_objects)} '{target_category}' object(s)...")

    for idx, obj in enumerate(filtered_objects):
        category = obj['category']
        bbox = obj['bbox']  # [x_min, y_min, x_max, y_max]

        # Extract bounding box coordinates
        x_min, y_min, x_max, y_max = bbox

        # Calculate bounding box dimensions
        bbox_width = x_max - x_min
        bbox_height = y_max - y_min

        # Make the bounding box square by expanding to the larger dimension
        max_dim = max(bbox_width, bbox_height)

        # Center the square around the original bbox center
        center_x = (x_min + x_max) / 2
        center_y = (y_min + y_max) / 2

        # Calculate square bbox
        x_min_square = center_x - max_dim / 2
        y_min_square = center_y - max_dim / 2
        x_max_square = center_x + max_dim / 2
        y_max_square = center_y + max_dim / 2

        # Now apply 25% padding to the square bbox
        pad = max_dim * (padding_percent / 100)

        # Apply padding and clamp to image boundaries
        x_min_padded = max(0, x_min_square - pad)
        y_min_padded = max(0, y_min_square - pad)
        x_max_padded = min(img_width, x_max_square + pad)
        y_max_padded = min(img_height, y_max_square + pad)

        # Crop the image
        cropped_img = img.crop((x_min_padded, y_min_padded, x_max_padded, y_max_padded))

        # Add scale bar (10 μm, 8 pixels/μm = 80 pixels)
        cropped_img = add_scale_bar(cropped_img, scale_bar_um=10, pixels_per_um=8)

        # Create output filename
        # Format: original_filename.png (no category suffix since we only crop stoma)
        output_filename = f"{base_name}.png"
        output_path = os.path.join(output_dir, output_filename)

        # Save cropped image with maximum quality and preserve DPI
        dpi = img.info.get('dpi', (96, 96))
        cropped_img.save(output_path, dpi=dpi, compress_level=0)
        print(f"  Saved: {output_filename} (size: {cropped_img.size}, DPI: {dpi})")


def crop_to_square_center(image_path, output_dir, scale_bar_um=10, pixels_per_um=11):
    """
    Crop an image to a square based on the shortest edge (centered).
    Then add a scale bar.

    Args:
        image_path: Path to the input image
        output_dir: Directory to save cropped image
        scale_bar_um: Length of scale bar in micrometers (default: 10)
        pixels_per_um: Pixels per micrometer (default: 11)
    """
    # Load the image
    img = Image.open(image_path)
    img_width, img_height = img.size

    print(f"Processing {image_path.name} (size: {img_width}x{img_height})...")

    # Find the shortest edge
    min_dim = min(img_width, img_height)

    # Calculate center
    center_x = img_width / 2
    center_y = img_height / 2

    # Calculate crop box for centered square
    x_min = center_x - min_dim / 2
    y_min = center_y - min_dim / 2
    x_max = center_x + min_dim / 2
    y_max = center_y + min_dim / 2

    # Crop to square
    cropped_img = img.crop((x_min, y_min, x_max, y_max))

    # Add scale bar
    cropped_img = add_scale_bar(cropped_img, scale_bar_um=scale_bar_um, pixels_per_um=pixels_per_um)

    # Create output filename (preserve original name but change extension to PNG)
    output_filename = image_path.stem + ".png"
    output_path = output_dir / output_filename

    # Save cropped image with maximum quality
    dpi = img.info.get('dpi', (96, 96))
    cropped_img.save(output_path, dpi=dpi, compress_level=0)
    print(f"  Saved: {output_filename} (size: {cropped_img.size}, DPI: {dpi})")


def process_stomata_annotation_set(input_dir, output_dir, padding_percent=25, target_category="stoma"):
    """
    Process a folder of ISAT JSON + PNG screenshot pairs and save cropped objects.

    Args:
        input_dir: Directory containing annotation JSON files and source PNG images
        output_dir: Directory to save cropped image outputs
        padding_percent: Percentage of padding around square bbox (default: 25)
        target_category: Object category to crop (default: "stoma")
    """
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Input directory: {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Padding: {padding_percent}%\n")

    json_files = sorted(input_dir.glob("*.json"))

    if not json_files:
        print("No JSON annotation files found!")
        return

    print(f"Found {len(json_files)} annotation file(s)\n")

    for json_path in json_files:
        try:
            crop_objects_from_annotation(
                json_path,
                output_dir,
                padding_percent=padding_percent,
                target_category=target_category
            )
        except Exception as e:
            print(f"Error processing {json_path.name}: {e}")

    print(f"\nDone! Cropped images saved to: {output_dir}")


def main():
    # Define paths
    base_dir = Path(__file__).parent

    # Process stomata aperture images
    screenshots_dir = base_dir / "data" / "Stomata aperture" / "examples" / "screen_shots"
    stomata_output_dir = base_dir / "data" / "Stomata aperture" / "examples"

    # Create output directory if it doesn't exist
    stomata_output_dir.mkdir(parents=True, exist_ok=True)

    print("=" * 80)
    print("PROCESSING STOMATA APERTURE IMAGES")
    print("=" * 80)
    process_stomata_annotation_set(
        screenshots_dir,
        stomata_output_dir,
        padding_percent=25,
        target_category="stoma"
    )

    # Process stomata width screenshots for supplementary figure
    print("\n" + "=" * 80)
    print("PROCESSING STOMATA WIDTH SCREENSHOTS")
    print("=" * 80)

    width_screenshots_dir = base_dir / "data" / "Stomata aperture" / "examples" / "width_sreen_shots"
    width_output_dir = width_screenshots_dir / "cropped"

    process_stomata_annotation_set(
        width_screenshots_dir,
        width_output_dir,
        padding_percent=25,
        target_category="stoma"
    )

    # Process blended peels starch staining images
    print("\n" + "=" * 80)
    print("PROCESSING BLENDED PEELS STARCH STAINING IMAGES")
    print("=" * 80)

    blended_input_dir = base_dir / "data" / "BL blended peels starch staining" / "examples" / "raw"
    blended_output_dir = base_dir / "data" / "BL blended peels starch staining" / "examples"

    # Create output directory if it doesn't exist
    blended_output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Input directory: {blended_input_dir}")
    print(f"Output directory: {blended_output_dir}")
    print("Scale bar: 10 um at 11 pixels/um\n")

    # Find all TIF files
    tif_files = list(blended_input_dir.glob("*.tif")) + list(blended_input_dir.glob("*.tiff"))

    if not tif_files:
        print("No TIF image files found!")
    else:
        print(f"Found {len(tif_files)} image file(s)\n")

        # Process each TIF file
        for tif_path in tif_files:
            try:
                crop_to_square_center(tif_path, blended_output_dir, scale_bar_um=10, pixels_per_um=11)
            except Exception as e:
                print(f"Error processing {tif_path.name}: {e}")

        print(f"\nDone! Cropped blended peels images saved to: {blended_output_dir}")

    # Process live/dead staining example images
    print("\n" + "=" * 80)
    print("PROCESSING LIVE/DEAD STAINING IMAGES")
    print("=" * 80)

    live_dead_input_dir = base_dir / "data" / "Live dead staining" / "examples" / "raw"
    live_dead_output_dir = base_dir / "data" / "Live dead staining" / "examples"

    live_dead_output_dir.mkdir(parents=True, exist_ok=True)

    print(f"Input directory: {live_dead_input_dir}")
    print(f"Output directory: {live_dead_output_dir}")
    print("Scale bar: 50 um at 2.8 pixels/um\n")

    live_dead_files = []
    for ext in ("*.png", "*.jpg", "*.jpeg", "*.tif", "*.tiff"):
        live_dead_files.extend(live_dead_input_dir.glob(ext))

    if not live_dead_files:
        print("No image files found!")
    else:
        print(f"Found {len(live_dead_files)} image file(s)\n")

        for image_path in live_dead_files:
            try:
                crop_to_square_center(image_path, live_dead_output_dir,
                                      scale_bar_um=50, pixels_per_um=2.8)
            except Exception as e:
                print(f"Error processing {image_path.name}: {e}")

        print(f"\nDone! Cropped live/dead images saved to: {live_dead_output_dir}")

    print("\n" + "=" * 80)
    print("ALL PROCESSING COMPLETE")
    print("=" * 80)


if __name__ == "__main__":
    main()
