ngenes=250

model=densenet121

# window=10
# window=30
window=74
# window=75
# window=224

# epochs=5
# epochs=20
# epochs=50
epochs=100

folds=3
# folds=4

# GPU=0
GPU=1
# GPU=2

# for patient in `python3 -m stnet patients`
# do
#      echo ${patient}
# done

# for patient in "BC30001"
# for patient in "BC30002"
# for patient in "BC30003"
# for patient in "BC30004"
# for patient in "BC30005"
for patient in "BC50027"
# for patient in "BC50027" "BC50040" "BC50111" "BC51218" "BC51517" "BC52337" "BC53934"
do
     echo ${patient}
     echo "Model : ${model}"
     echo "Epochs: ${epochs}"
     echo "Window: ${window}px"
     echo "GPU   : ${GPU}"

     CUDA_VISIBLE_DEVICES=${GPU} bin/cross_validate.py output/${model}_${window}/top_${ngenes}/${patient}_ ${folds} ${epochs} ${patient} --lr 1e-6 --window ${window} --model ${model} --pretrain --average --batch 32 --workers 7 --gene_n ${ngenes} --norm --epochs ${epochs} 
done