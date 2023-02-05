
# 开启Python虚拟环境
source activate py376

# 定位到输入文件所在位置
cd /Users/mosquito/Desktop/NAD+\ NetBID/NAD/project_2022-09-21/SJAR/project_2022-09-21/

# TF
sjaracne local -e input.exp -g tf.txt -n 100 -o ./TF_out -pc 1e-5

# SIG
sjaracne local -e input.exp -g sig.txt -n 100 -o ./SIG_out -pc 1e-5

