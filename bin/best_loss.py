#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK

import argparse
import argcomplete

parser = argparse.ArgumentParser()
parser.add_argument("--epochs", type=int)
parser.add_argument("--npz_path", type=str)

argcomplete.autocomplete(parser)
args = parser.parse_args()

import os
# import sys
from PIL import Image
from glob import glob
from natsort import natsorted
import pathlib
import numpy as np
import matplotlib.pyplot as plt

if __name__ == "__main__":

    # args.npz_path = "/home/tester/bokleong/ST-Net/output/train_COVIDHCC_test_BC30004"

    npzs = natsorted(glob(f"{args.npz_path}/*.npz", recursive=True))

    losses = []
    best_loss = float("inf")

    for epoch, npz in enumerate(npzs):

        if not "_cv.npz" in npz:
            data = np.load(npz)

            patient = []
            pred = []
            count = []
            mean = []

            pred.append(data["predictions"])
            count.append(data["counts"])
            mean.append(data["mean_expression"])
            patient.append(data["patient"])

            mean = np.concatenate([np.repeat(np.expand_dims(m, 0), c.shape[0], 0) for (m, c) in zip(mean, count)], 0)
            pred = np.concatenate(pred)
            count = np.concatenate(count)
            patient = np.concatenate(patient)

            loss = np.sum((pred - count) ** 2)
            losses.append(loss)

            if loss < best_loss:
                best_loss = loss
                best_epoch = epoch + 1

    fig, ax = plt.subplots(figsize=(20,8))
    ax.plot(losses)
    ax.set_xticks(np.arange(1, 101))
    ax.set_xticklabels(np.arange(1, 101), rotation=45, ha="right")
    ax.grid(which="major", axis="x")
    # ax.minorticks_on()
    # ax.grid(which="minor", axis="x")
    plt.tight_layout()
    plt.savefig(f"{args.npz_path}/00-losses.png")