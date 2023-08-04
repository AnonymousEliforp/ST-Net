#!/usr/bin/env python3
"""
Calculate Spearman's correlation for ST-Net predictions
"""
import os
import argparse
import argcomplete
parser = argparse.ArgumentParser()
parser.add_argument("output_dir", type=str, nargs="+")
parser.add_argument("--epoch", type=str, default="cv")
parser.add_argument("--plotname", type=str, default="boxplot")
parser.add_argument("--plotlims", type=float, nargs=2, default=[0.0,0.0])
argcomplete.autocomplete(parser)
args = parser.parse_args()
import sys
sys.path.append(".")
import json
import stnet
import pathlib
import scipy.stats
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from matplotlib.ticker import MultipleLocator

if __name__ == "__main__":

    # print(f"{args.epoch = }")
    npzs = glob(f"{args.output_dir[0]}/*_{args.epoch}.npz",recursive=True)
    # Sort npzs found by glob, so that boxplot is logical and pleasant
    npzs.sort(key=lambda x: int(os.path.basename(x).split(f"_{args.epoch}.")[0][2:]))
    # print(f"{args.output_dir[0] = }")
    # print(f"{npzs = }")
    # print(f"{len(npzs) = }")

    patient_map = {
        "BC30001":"COVID_HCC1",
        "BC30002":"COVID_HCC2",
        "BC30003":"COVID_HCC3",
        "BC30004":"COVID_HCC4",
        "BC30005":"COVID_HCC5",
        "BC30006":"COVID_HCC6",
        "BC30007":"COVID_HCC7",
        "BC42334":"NCC023Pre_NCC034Pre",
        "BC42740":"NCC027Pre_NCC040A2Pre",
        "BC43740":"NCC037Pre_NCC040A1Pre",
        "BC50027":"NCC027Post",
        "BC50040":"NCC040Post",
        "BC50111":"NCC001Post_NCC011Post",
        "BC51218":"NCC012Post_NCC018Post",
        "BC51517":"NCC015Post_NCC017Post",
        "BC52337":"NCC023Post_NCC037Post",
        "BC53934":"NCC039Post_NCC034Post",
        "HCC12"  :"COVID_HCC1,2",
        "HCC1234":"COVID_HCC1,2,3,4"
    }

    vals = {}
    all_spearmanr = {}
    for npz in npzs:
        
        data = np.load(npz,allow_pickle=True)
        patient     = patient_map[data["patient"][0]]
        counts      = data["counts"]
        predictions = data["predictions"]
        gene_names  = data["gene_names"]
        this_npz_vals = []
        save_spearmanr = {}

        for gene in gene_names:
            index = np.argmax(gene_names == gene)
            # print(f"{index = }")
            c = counts[:,index]
            p = predictions[:,index]
            spearmanr, pval = scipy.stats.spearmanr(c,p)
            
            # Only append 5% significant correlations
            if pval < 0.05:
                this_npz_vals.append(spearmanr)
                save_spearmanr[gene] = spearmanr

        vals[patient] = this_npz_vals
        save_spearmanr = {k: v for k, v in sorted(save_spearmanr.items(), key=lambda item: item[1], reverse=True)}
        all_spearmanr[patient] = save_spearmanr

    # print(vals)
    # print(len(vals))
    # print(len(vals[0]))
    # print(len(vals[1]))
    print(f"{gene_names=}")
    print(f"{predictions=}")
    print(f"{counts=}")
    

    title_size = 20
    label_size = 16
    fig,ax = plt.subplots(1,1,figsize=(12,8))
    # print(f"{vals = }")
    ax.boxplot(vals.values(),
               labels=vals.keys(),
               patch_artist=True)
    model_testpatient = patient_map[args.plotname.split('_')[0]]
    ax.set_title(f"Spearman's Correlations for inference using\n{model_testpatient}-trained model\n",fontsize=title_size)
    ax.set_xlabel(f"\nPatients",fontsize=label_size)
    ax.set_xticklabels(vals.keys(),rotation=45,ha='right')
    ax.set_ylabel(f"Spearman's Correlation\n",fontsize=label_size)
    ax.tick_params(axis='both', which='major', labelsize=label_size)
    ax.minorticks_on()
    ax.grid(which="minor",axis='y',linewidth=0.3,linestyle='--')
    ax.grid(which="major",axis='y',linewidth=1.2)
    ax.yaxis.set_major_locator(MultipleLocator(0.1))
    ax.axhline(y=0,color='r')
    
    if not args.plotlims == [0.0,0.0]:
        print(f"{args.plotlims = }")
        lower = args.plotlims[0]
        upper = args.plotlims[1]
        ax.set_ybound(lower=lower,upper=upper)

    # plt.tight_layout()

    # Add scatter
    npatients = len(vals)
    print(f"vals=\n{vals}")
    """
    for i, d in enumerate(data):
    """
    ALB_CD74_coord_list = [] 
    vals_new = list(list(vals.values())[0])
    print(f"{vals_new=}")
    top_i = sorted(range(len(vals_new)), key=lambda i: vals_new[i], reverse=True)[:10]
    print(vals_new)
    direction = 1
    y_flag = 0.0
    for i in top_i:
        #if name=="ALB" or name=="CD74":
        ALB_CD74_coord_list.append((gene_names[i],(1,i)))

        # Coordinates of the point
        point_x = 1.0
        point_y = list(list(vals.values())[0])[i]
        print("point_y", point_y)
        # Add the label with an arrow and text
        plt.annotate(gene_names[i], xy=(point_x, point_y), xytext=(point_x + 0.2*(direction), point_y - y_flag), arrowprops=dict(arrowstyle='->'), fontsize=30)
        direction *= -1
        y_flag += 0.05


    plt.scatter([1]* len(list(vals.values())[0]),list(vals.values())[0], c = 'red', zorder = 100)
    plt.show()
    
    # figroot = args.output_dir[0].split('/')[:-3]
    # figroot = '/'.join(figroot)
    # figroot = f"{figroot}/fig/visualize/"
    pathlib.Path(os.path.dirname(args.output_dir[0])).mkdir(parents=True, exist_ok=True)
    plt.savefig(f"{args.output_dir[0]}/{args.plotname}.png",dpi=300, bbox_inches="tight")
    print(f"Saved to {args.output_dir[0]}/{args.plotname}.png")

    with open(f"{args.output_dir[0]}/{args.plotname}.json", "w") as f:
        json.dump(save_spearmanr, f, indent=4)
