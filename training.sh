ngenes=250

model=densenet121

# window=10
window=296
# window=74
# window=75
# window=224

# epochs=5
# epochs=20
# epochs=50
epochs=100

folds=3
# folds=4

GPU=0
# GPU=1
# GPU=2

using_to_train="COVIDHCC"
# using_to_train="Y90Pre"
# using_to_train="Y90Post"

# for patient in "BC30001"
# for patient in "BC30002"
# for patient in "BC30003"
#change this based on the patients that you want to eval on 
for patient in "BC30007"
# for patient in "BC30005"
# for patient in "BC50027"
# for patient in "BC50111"
# for patient in "BC50027" "BC50040" "BC50111" "BC51218" "BC51517" "BC52337" "BC53934"
do
     echo ${patient}
     echo "Model : ${model}"
     echo "Epochs: ${epochs}"
     echo "Window: ${window}px"
     echo "GPU   : ${GPU}"
     echo "Using : ${using_to_train}"

     CUDA_VISIBLE_DEVICES=${GPU} bin/cross_validate.py output/train_${using_to_train}_test_${patient}/${patient}_ ${folds} ${epochs} ${patient} --lr 1e-6 --window ${window} --model ${model} --pretrain --average --batch 32 --workers 7 --gene_n ${ngenes} --norm --epochs ${epochs} 
     
done
