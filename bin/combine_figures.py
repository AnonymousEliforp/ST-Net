#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

import argparse
import argcomplete

parser = argparse.ArgumentParser()
parser.add_argument("--patient", type=str)
parser.add_argument("--gene", type=str)
parser.add_argument("--model", type=str)
parser.add_argument("--inference_type", type=str)
parser.add_argument("--figroot", default="fig/", type=str)

argcomplete.autocomplete(parser)
args = parser.parse_args()

import os
# import sys
from PIL import Image
import pathlib
pathlib.Path(os.path.dirname(args.figroot)).mkdir(parents=True, exist_ok=True)

gene_order_dict = {
        "ALB": 1,
        "SERPINA1": 2,
        "APOA1": 3,
        "APOC1": 4,
        "IGKC": 5,
        "APOA2": 6,
        "FGB": 7,
        "FGA": 8,
        "C1S": 9,
        "CD74": 10,
        "MT2A": 11
    }

pred_image = f"output/train_{args.inference_type}_test_{args.patient}/{args.patient}_visualize/[{gene_order_dict[args.gene]}] {args.patient}_{args.gene}_pred.png"
gt_image   = f"output/train_{args.inference_type}_test_{args.patient}/{args.patient}_visualize/[{gene_order_dict[args.gene]}] {args.patient}_{args.gene}.png"
scatter_image   = f"output/train_{args.inference_type}_test_{args.patient}/{args.patient}_visualize/[{gene_order_dict[args.gene]}] {args.patient}_{args.gene}_scatter.png"
ranked_image   = f"output/train_{args.inference_type}_test_{args.patient}/{args.patient}_visualize/[{gene_order_dict[args.gene]}] {args.patient}_{args.gene}_ranked.png"

# print(f"{pred_image = }")
# print(f"{gt_image   = }")
# print(f"{scatter_image   = }")

# For COVID_HCC
crop = [
  (20, 30),
  (30, 20),
  (20, 50),
  (10, 50)
]

# For Y90Post
# crop = [
#   (60, 60),
#   (60, 60),
#   (80, 160),
#   (100, 160)
# ]

images = [Image.open(x) for x in [pred_image, gt_image, scatter_image, ranked_image]]
# print(f"{images = }")
cropped_images = []
for i, image in enumerate(images):
  width, height = image.size
  crop_dimensions = (crop[i][0], 0, width - crop[i][0], height)
  cropped_images.append(image.crop(crop_dimensions))
# widths, heights = zip(*(i.size for i in images))
widths, heights = zip(*(i.size for i in cropped_images))

total_width = sum(widths)
max_height = max(heights)

new_im = Image.new('RGB', (total_width, max_height))

x_offset = 0
# for im in images:
for im in cropped_images:
  new_im.paste(im, (x_offset,0))
  x_offset += im.size[0]

pathlib.Path(os.path.dirname(f"{args.figroot}{args.patient}_{args.model}-model/")).mkdir(parents=True, exist_ok=True)
new_im.save(f"{args.figroot}{args.patient}_{args.model}-model/{args.patient}_{args.model}-model [{gene_order_dict[args.gene]}] {args.gene}_combined.png")