import os
import numpy as np
import pandas as pd
import pickle
import logging
import pathlib
import datetime
import time
import glob
import collections
import stnet
import tqdm


def spatial(args):
    import skimage.io

    window = 224  # only to check if patch is off of boundary

    print(f"{args.logfile = }")
    print(f"{args.loglevel = }")
    stnet.utils.logging.setup_logging(args.logfile, args.loglevel)

    logger = logging.getLogger(__name__)

    pathlib.Path(args.dest).mkdir(parents=True, exist_ok=True)

    raw, subtype = load_raw(args.root)

    logger.info(f"{len(raw) = }")
    logger.info(f"{args.root = }")
    logger.info(f"{args.dest = }")
    logger.info(f"{subtype = }")

    with open(args.dest + "/subtype.pkl", "wb") as f:
        pickle.dump(subtype, f)

    t = time.time()

    t0 = time.time()
    section_header = None
    gene_names = set()
    for patient in raw:
        for section in raw[patient]:
            section_header = raw[patient][section]["count"].columns.values[0]
            gene_names = gene_names.union(set(raw[patient][section]["count"].columns.values[1:]))
    gene_names = list(gene_names)
    gene_names.sort()
    print(f"{section_header = }")
    # print(f"{gene_names = }")
    with open(args.dest + "/gene.pkl", "wb") as f:
        pickle.dump(gene_names, f)
    gene_names = [section_header] + gene_names
    # print(f"{gene_names = }")
    logger.info("Finding list of genes: " + str(time.time() - t0))

    for (i, patient) in enumerate(raw):
        logger.info("Processing " + str(i + 1) + " / " + str(len(raw)) + ": " + patient)

        for section in raw[patient]:

            pathlib.Path("{}/{}/{}".format(args.dest, subtype[patient], patient)).mkdir(parents=True, exist_ok=True)

            # This is just a blank file to indicate that the section has been completely processed.
            # Preprocessing occassionally crashes, and this lets the preparation restart from where it let off
            complete_filename = "{}/{}/{}/.{}".format(args.dest, subtype[patient], patient, section)
            if pathlib.Path(complete_filename).exists():
                logger.info("Patient {} section {} has already been processed.".format(patient, section))
            else:
                logger.info("Processing " + patient + " " + section + "...")

                # In the original data, genes with no expression in a section are dropped from the table.
                # This adds the columns back in so that comparisons across the sections can be done.
                t0 = time.time()
                missing = list(set(gene_names) - set(raw[patient][section]["count"].keys()))
                # print(f"{missing = }")
                print(f"{len(missing) = }")
                c = raw[patient][section]["count"].values[:, 1:].astype(float)
                # Adding columns for dropped genes, and setting all values in that column to zero
                pad = np.zeros((c.shape[0], len(missing)))
                c = np.concatenate((c, pad), axis=1)
                names = np.concatenate((raw[patient][section]["count"].keys().values[1:], np.array(missing)))
                # Ordering columns?
                c = c[:, np.argsort(names)]
                logger.info("Adding zeros and ordering columns: " + str(time.time() - t0))

                t0 = time.time()
                count = {}
                for (index, row) in raw[patient][section]["count"].iterrows():
                    count[row.values[0]] = c[index, :]
                logger.info("Extracting counts: " + str(time.time() - t0))

                t0 = time.time()
                tumor = {}
                not_int = False
                for (_, row) in raw[patient][section]["tumor"].iterrows():
                    if isinstance(row[1], float) or isinstance(row[2], float):
                        not_int = True
                    # print(f"{row.tolist() = }")
                    x = int(round(row[1]))
                    y = int(round(row[2]))
                    # logger.info(f"tumor: x = {x:<3}, y = {y:<3}")
                    tumor[(x, y)] = (row[4] == "tumor")
                if not_int:
                    logger.warning("Patient " + patient + " " + section + " has non-integer patch coordinates.")
                logger.info("Extracting tumors: " + str(time.time() - t0))
                print(raw[patient][section]["image"])
                t0 = time.time()
                image = skimage.io.imread(raw[patient][section]["image"])
                logger.info("Loading image: " + str(time.time() - t0))

                data = []
                for (_, row) in raw[patient][section]["spot"].iterrows():
                    # x = int(round(row["pixel_x"]))
                    # y = int(round(row["pixel_y"]))
                    # x = int(round(row["X"]))
                    # y = int(round(row["Y"]))

                    # print(f"{row = }")
                    # print(f"{type(row) = }")
                    # print(f"{row.tolist()[0].split(',') = }")
                    tmp_list = row.tolist()[0].split(',')
                    pixel_x = int(round(float(tmp_list[1])))
                    pixel_y = int(round(float(tmp_list[2])))

                    x_coord = int(tmp_list[0].split('x')[0])
                    y_coord = int(tmp_list[0].split('x')[1])
                    # print(f"{x_coord = }")
                    # print(f"{y_coord = }")

                    X = image[(pixel_y + (-window // 2)):(pixel_y + (window // 2)), (pixel_x + (-window // 2)):(pixel_x + (window // 2)), :]
                    # if X.shape == (window, window, 3):
                    if True:
                        # if (int(row["X"]), int(row["Y"])) in tumor:
                        if (x_coord, y_coord) in tumor:
                            
                            vals = np.unique(count[str(x_coord) + "x" + str(y_coord)])
                            num_nan = np.sum(np.isnan(vals))

                            if num_nan == 0:
                                logger.info(f"Appending {x_coord}x{y_coord}...")
                                data.append((X,
                                            #  count[str(int(row["x"])) + "x" + str(int(row["y"]))],
                                            count[str(x_coord) + "x" + str(y_coord)],
                                            #  tumor[(int(row["x"]), int(row["y"]))],
                                            tumor[(x_coord, y_coord)],
                                            #  np.array([x, y]),
                                            np.array([pixel_x, pixel_y]),
                                            np.array([patient]),
                                            np.array([section]),
                                            #  np.array([int(row["x"]), int(row["y"])]),
                                            np.array([x_coord, y_coord])
                                            ))

                                # filename = "{}/{}/{}/{}_{}_{}.npz".format(args.dest, subtype[patient], patient, section,
                                #                                           int(row["x"]), int(row["y"]))
                                filename = "{}/{}/{}/{}_{}_{}.npz".format(args.dest, subtype[patient], patient, section,
                                                                        x_coord, y_coord)
                                
                                # np.savez_compressed(filename, count=count[str(int(row["x"])) + "x" + str(int(row["y"]))],
                                #                     tumor=tumor[(int(row["x"]), int(row["y"]))],
                                #                     pixel=np.array([x, y]),
                                #                     patient=np.array([patient]),
                                #                     section=np.array([section]),
                                #                     index=np.array([int(row["x"]), int(row["y"])]))
                                np.savez_compressed(filename, count=count[str(x_coord) + "x" + str(y_coord)],
                                                    tumor=tumor[(x_coord, y_coord)],
                                                    pixel=np.array([pixel_x, pixel_y]),
                                                    patient=np.array([patient]),
                                                    section=np.array([section]),
                                                    index=np.array([x_coord, y_coord]))

                            else:
                                logger.info(f"NaN values in count. Not appending {x_coord}x{y_coord}.")

                        else:
                            logger.warning("Patch " + str(x_coord) + "x" + str(
                                y_coord) + " not found in " + patient + " " + section)
                    else:
                        logger.warning("Detected spot too close to edge.")
                        pass
                logger.info("Saving patches: " + str(time.time() - t0))

                with open(complete_filename, "w"):
                    pass
    logger.info("Preprocessing took " + str(time.time() - t) + " seconds")

    if (not os.path.isfile(stnet.config.SPATIAL_PROCESSED_ROOT + "/mean_expression.npy") or
        not os.path.isfile(stnet.config.SPATIAL_PROCESSED_ROOT + "/median_expression.npy")):
        logging.info("Computing statistics of dataset")
        gene = []
        logger.info(f"{args.dest = }")
        for filename in tqdm.tqdm(glob.glob("{}/*/*/*_*_*.npz".format(args.dest))):
            npz = np.load(filename)
            count = npz["count"]
            gene.append(np.expand_dims(count, 1))

        gene = np.concatenate(gene, 1)
        np.save(stnet.config.SPATIAL_PROCESSED_ROOT + "/mean_expression.npy", np.mean(gene, 1))
        np.save(stnet.config.SPATIAL_PROCESSED_ROOT + "/median_expression.npy", np.median(gene, 1))


def load_section(root: str, patient: str, section: str, subtype: str):
    """
    Loads data for one section of a patient.
    """
    import pandas
    import gzip

    # file_root = root + "/" + subtype + "/" + patient + "/" + patient + "_" + section
    # file_root = root + "/" + "BC" + patient + "_" + section
    file_root = root
    BC = "BC" + patient[2:] + "_" + section
    BT = patient + "_" + section


    # image = skimage.io.imread(file_root + ".jpg")
    count_matrix = file_root + BC + "_stdata"
    # image = file_root + "HE_" + BT + ".jpg"
    image = file_root + "HE_" + BT + ".tif"
    spot_coordinates = file_root + "spots_" + BT
    tumor_annotation = file_root + BC + "_Coords"

    # print(f"\n{file_root = }")
    # print(f"{count_matrix = }")
    # print(f"{image = }")
    # print(f"{spot_coordinates = }")
    # print(f"{tumor_annotation = }")

    # Count matrix
    # if stnet.utils.util.newer_than(file_root + ".tsv.gz", file_root + ".pkl"):
    if stnet.utils.util.newer_than(count_matrix + ".tsv.gz", count_matrix + ".pkl"):
        # with gzip.open(file_root + ".tsv.gz", "rb") as f:
        with gzip.open(count_matrix + ".tsv.gz", "rb") as f:
            count = pandas.read_csv(f, sep="\t")
        with open(count_matrix + ".pkl", "wb") as f:
            pickle.dump(count, f)
    else:
        with open(count_matrix + ".pkl", "rb") as f:
            count = pickle.load(f)

    # Spot coordinates
    # if stnet.utils.util.newer_than(file_root + ".spots.txt", file_root + ".spots.pkl"):
    if stnet.utils.util.newer_than(spot_coordinates + ".csv.gz", spot_coordinates + ".pkl"):
        # spot = pandas.read_csv(file_root + ".spots.txt", sep="\t")
        with gzip.open(spot_coordinates + ".csv.gz", "rb") as f:
            spot = pandas.read_csv(f, sep="\t")
        with open(spot_coordinates + ".pkl", "wb") as f:
            pickle.dump(spot, f)
    else:
        with open(spot_coordinates + ".pkl", "rb") as f:
            spot = pickle.load(f)

    # Tumor annotation
    if stnet.utils.util.newer_than(tumor_annotation + ".tsv.gz", tumor_annotation + ".pkl"):
        # tumor = pandas.read_csv(tumor_annotation + "_Coords.tsv", sep="\t")
        with gzip.open(tumor_annotation + ".tsv.gz", "rb") as f:
            tumor = pandas.read_csv(f, sep="\t")
        with open(tumor_annotation + ".pkl", "wb") as f:
            pickle.dump(tumor, f)
    else:
        with open(tumor_annotation + ".pkl", "rb") as f:
            tumor = pickle.load(f)

    return {"image": image, "count": count, "spot": spot, "tumor": tumor}


def load_raw(root: str):
    """
    Loads data for all patients.
    """
    # __location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__))) + "/"
    print(f"{root = }")

    logger = logging.getLogger(__name__)

    # Wildcard search for patients/sections
    # images = glob.glob(__location__ + root + "*_*_*.jpg")
    # images = glob.glob(root + "*_*.jpg")
    images = glob.glob(root + "*_*.tif")
    df = pd.read_csv(root + "metadata.csv")
    print(f"\ndf.head():\n{df.head()}\n")

    # print(f"{type(images) = }")
    # print(f"{len(images) = }")
    # print(f"{images[-1] = }")

    temp = []
    for image in images:
        this_string = image.replace("\\", "/")
        temp.append(this_string)

    images = temp
    print(images)

    # Dict mapping patient ID (str) to a list of all sections available for the patient (List[str])
    patient = collections.defaultdict(list)
    for (useless, p, s) in map(lambda x: x.split("/")[-1][:-4].split("_"), images):
        # print(f"{p = }")
        # print(f"{s = }")
        patient[p].append(s)

    # Dict mapping patient ID (str) to subtype (str)
    subtype = {}
    # for (st, p) in map(lambda x: (x.split("/")[-3], x.split("/")[-1][:-4].split("_")[1]), images):
    for p in map(lambda x: x.split("/")[-1][:-4].split("_")[1], images):
        # print(f"{p = }")
        p2 = "BC" + p[2:]
        print( p2, df["patient"])
        st = df.loc[df["patient"] == p2, "type"].values[0]
        # print(f"{type(st) = }")
        # print(f"st = {st}")
        if p in subtype:
            # print(f"Already in subtype\n")
            # print(f"[Y]\n")
            if subtype[p] != st:
                raise ValueError("Patient {} is marked as type {} and {}.".format(p, subtype[p], st))
        else:
            # print(f"Not in subtype\n")
            # print(f"[N]\n")
            subtype[p] = st

    logger.info("Loading raw data...")
    t = time.time()
    data = {}
    with tqdm.tqdm(total=sum(map(len, patient.values()))) as pbar:
        for p in patient:
            data[p] = {}
            for s in patient[p]:
                data[p][s] = load_section(root, p, s, subtype[p])
                pbar.update()
    logger.info("Loading raw data took " + str(time.time() - t) + " seconds.")

    # print(f"{data = }")

    return data, subtype
