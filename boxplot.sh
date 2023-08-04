#!/bin/bash

## [1] SET THIS
## Select the inference patient type
inference_type="COVIDHCC"
# inference_type="Y90Pre"
# inference_type="Y90Post"

## [2] SET THIS
## Select the model used for inference
# model=BC30001
model=BC30007
# model=BC50027
# model=BC50040
# model=BC50111
# model=BC51218
# model=BC51517
# model=BC52337
# model=BC53934
# model=HCC12
# model=HCC1234

## [3] SET THIS
window=30

# for model in "BC50027" "BC50040" "BC50111" "BC51218" "BC51517" "BC52337" "BC53934"
for model in "BC30007"
do
    ## Set automatically
    output_dir="${PWD}/output/train_${inference_type}_test_${model}"
    # output_dir="$PWD/output/densenet121_$window/top_250"
    epoch=cv
    # epoch=50
    # epoch=100

    # bin/plot_boxplots.py ${output_dir} --epoch ${epoch} --plotname ${model}_model_inference --plotlims -0.5 0.5
    bin/plot_boxplots.py ${output_dir} --epoch ${epoch} --plotname ${model}_model_inference_epoch_${epoch}
done
