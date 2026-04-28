import os
from PIL import Image
from reportlab.pdfgen import canvas

#process png to pdf
def images_to_pdf(folder_path, output_pdf, include_title=False):
    # Collect all PNG files and sort by name
    image_files = sorted([f for f in os.listdir(folder_path) if f.endswith('.png')])

    # Create a new PDF file
    c = canvas.Canvas(output_pdf)

    for image_file in image_files:
        img_path = os.path.join(folder_path, image_file)
        img = Image.open(img_path)

        # Retrieve image dimensions and convert to PDF units (pt)
        width, height = img.size
        width, height = width * 0.75, height * 0.75  # assuming 1 pixel = 0.75 pt

        # Set the page size
        if include_title:
            c.setPageSize((width, height + 50))  # Reserve space for the title
            c.drawString(10, height + 20, image_file)  # Draw the title at the top
        else:
            c.setPageSize((width, height))

        # Draw the image
        c.drawImage(img_path, 0, 0, width=width, height=height)
        c.showPage()  # Start a new page

    c.save()

# folder_path = "source_data/Robin/Robin_png/"  # Replace with your folder path
# output_pdf = "Robin_RNA.pdf"  # Path for the output PDF file
# include_title = True  # Whether to include a title on each page
# images_to_pdf(folder_path, output_pdf, include_title)


from PIL import Image
#need pip install pillow

def convert_white_to_transparent(input_image_path, output_image_path):
    """
    Transfer graph with white color as background color into opticify as background color.
    """
    image = Image.open(input_image_path).convert("RGBA")

    # Retrieve the image's pixel data
    datas = image.getdata()

    new_data = []
    for item in datas:
        # Replace white (and near-white) pixels with transparency
        if item[0] > 200 and item[1] > 200 and item[2] > 200:
            new_data.append((255, 255, 255, 0))  # Set pixel to fully transparent
        else:
            new_data.append(item)

    # Update the image's pixel data
    image.putdata(new_data)

    # Save the modified image
    image.save(output_image_path, "PNG")


"""
Example Usage:
input_image_path = "pic/Pic_new/SAM_II.png"  # input filepath of graph
output_image_path = "pic/Pic_new/SAM_II_new.png"  # output filepath of graph

convert_white_to_transparent(input_image_path, output_image_path)
print(f"Saved transparent image to {output_image_path}")
"""

