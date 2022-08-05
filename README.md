<!--
 * @Author: haoqiang haoqiang@mindrank.ai
 * @Date: 2022-07-19 03:41:19
 * @LastEditors: haoqiang haoqiang@mindrank.ai
 * @LastEditTime: 2022-08-05 08:56:25
 * @FilePath: /work-home/molecule-3d-similarity/README.md
 * @Description: 
 * 
 * Copyright (c) 2022 by haoqiang haoqiang@mindrank.ai, All Rights Reserved. 
-->
# 1 Overview
| 序号 | 名称 | 来源/下载地址 | 相关文献 | 安装步骤 | 测试 |
| -- | -- | -- | -- | -- | -- |
| 01 | ShapeProtrudeDist | rdkit | - | - | ✓ |
| 02 | ShapeTanimotoDist | rdkit | - | - | ✓ |
| 03 | ShapeTverskyIndex | rdkit | - | - | ✓ |
| 04 | SC score | rdkit | [@Yang2020](https://pubs.rsc.org/en/content/articlelanding/2020/sc/d0sc03126g) | - | ✓ |
| 05 | Gobbi_Pharm2D | rdkit | - | - | ✓ |
| 06 | USR | oddt | [@Ballester2007](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.20681) | - | ✓ |
| 07 | USRCAT | oddt | [@Schreyer2012](https://doi.org/10.1186/1758-2946-4-27) | - | ✓ |
| 08 | Electroshape | oddt | [@Armstrong2010](https://doi.org/10.1007/s10822-010-9374-0) | - | ✓ |
| 09 | ACPC | [riken](http://www2.riken.jp/zhangiru/software.html) | [@Berenger2014](https://doi.org/10.1186/1758-2946-6-23) | [Installation](#ACPC) | ✓ |
| 10 | ESP-Sim | [pip](https://pypi.org/project/espsim/) | [@Bolcalto2021](https://chemrxiv.org/engage/chemrxiv/article-details/6182a7f68ac7a22cf566624d) | - | ✓ |
| 11 | PAPER | [SimTK](https://simtk.org/projects/paper) | [@Haque2010](https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.21307) | [software](#PAPER) | ✗ |
| 12 | MolShaCS | [code.google](https://code.google.com/archive/p/molshacs/downloads) | [@VazdeLima2013](https://www.sciencedirect.com/science/article/abs/pii/S0223523412006824) | - | ✗ |
| 13 | 🌟SHAFTS | [lilab](http://lilab-ecust.cn/home/resource.html) | [@Liu2011](https://doi.org/10.1021/ci200060s) | [software](#SHAFTS) | ✗ |
| 14 | ShaEP | [mivainio](http://users.abo.fi/mivainio/shaep/index.php) | [@Vainio2009](https://doi.org/10.1021/ci800315d) | - | ✗ |
| 15 | SimG | [lilab](http://lilab.ecust.edu.cn/home/resource.html) | [@Cai2013](https://doi.org/10.1021/ci400139j) | [software](#SimG) | ✗ |
| 16 | EGNN | [git](https://gitlab.mr.ai/haoqiang/egnn) | [@Satorras2022](http://arxiv.org/abs/2102.09844) | [usage](#EGNN) | ✗ |
|  |  |  |  | - | ✗ |


> - 相关综述：
>   - [[@Bero2017](https://iopscience.iop.org/article/10.1088/1742-6596/892/1/012015)] Similarity Measure for Molecular Structure: A Brief Review
>   - [[@Kumar2018](https://www.frontiersin.org/articles/10.3389/fchem.2018.00315)] Advances in the Development of Shape Similarity Methods and Their Application in Drug Discovery
>   - [[@Jiang2021a](https://doi.org/10.1093/bib/bbab231)] A comprehensive comparative assessment of 3D molecular similarity tools in ligand-based virtual screening

# 2 Installation
## 2.1 <span id='ACPC'>ACPC</span>
---
1. Install [Opam](https://opam.ocaml.org/doc/Install.html)
    ```shell
    sudo apt install software-properties-common
    sudo add-apt-repository ppa:avsm/ppa # may get message that `ERROR: '~avsm' user or team does not exist`, please try again
    sudo apt update
    sudo apt install gnuplot-x11 autoconf opam make gcc patch
    ```
2. Initialization Opam
    ```shell
    rm -rf ~/.opam
    opam init --disable-sandboxin # for LXC container
    ```
3. Install ACPC
    ```shell
    opam install ACPC # try until prompt change from `Processing  1/1: [default: http]` to `Processing  1/1`:
    ```
## 2.2 <span id='PAPER'>PAPER</span>
---
已下载在molecule-3d-similarity/softwares/目录下

## 2.3 <span id='SHAFTS'>SHAFTS</span>
---
已下载在molecule-3d-similarity/softwares/目录下, 如果出现`./Cynthia: No such file or directory`错误，可以尝试安装`sudo apt-get install lib32stdc++6`（[参考](https://haow.ca/blog//2017/No-such-file/)）。

## 2.4 <span id='SimG'>SimG</span>
---
已下载在molecule-3d-similarity/softwares/目录下

## 2.5 <span id='EGNN'>EGNN</span>
---
### 2.5.1 文件树
```shell
.
├── ae_datasets
│   ├── dataloader.py
│   ├── d_creator.py
│   ├── d_selector.py
│   └── __init__.py
├── dude # need to notice
│   ├── data # dude preprocess data
│   │   ├── aa2ar.csv
│   │   ├── abl1.csv
│   │   ├── ace.csv
│   │   ├── aces.csv
...
│   │   └── xiap.csv
│   ├── data.py # load data utils
│   ├── __init__.py
│   ├── models.py # construct models
│   └── utils.py
├── eval.py
├── graph.py
├── LICENSE
├── losess.py
├── main_ae.py
├── main_dude.py # dude main
├── main_nbody.py
├── main_qm9.py
├── models
│   ├── ae.py
│   ├── egnn_clean
│   │   ├── egnn_clean.py
│   │   └── __init__.py
│   ├── egnn.png
│   ├── gcl.py
│   └──  __init__.py
├── n_body_system
├── qm9
├── README.md
└── utils.py
```
### 2.5.2 模型结构
```
bs -> batch size
nn -> node number
nt -> node type
cp -> charge power
pm -> postion number (len(x, y, z) = 3)
```
1. 模型构建参数  

    |   | 名称       | 含义           | 类型         | 取值     | 备注                                                |
    | - | ---------- | -------------- | ------------ | -------- | --------------------------------------------------- |
    | 1 | in_node_nf | 输入节点类型数 | int          | 11*(2+1) | nf=number of feature；<br/>in_node_nf=nt\*(cp+1); |
    | 2 | in_edge_nf | 输入边类型数   | int          | 4        | (single, double, triple, aromatic)                  |
    | 3 | hidden_nf  | embedding维度  | int          | 128      |                                                     |
    | 4 | device     | 模型载入设备   | torch.device | -        |                                                     |
    | 5 | n_layers   | E_GCL数量      | int          | 7        |                                                     |
    | 6 | attention  | 注意力机制     | bool         | -        |                                                     |
2. 模型输入参数  
   
    |  | 名称      | 含义            | 形状                    | 关联中间变量                                             | 备注                                                                                                                        |
    | - | --------- | --------------- | ----------------------- | -------------------------------------------------------- | --------------------------------------------------------------------------------------------------------------------------- |
    | 1 | h         | Nodes输入向量   | (2, bs\*nn, nt\*(cp+1)) | one_hot, chargs, charge_power,<br />charge_scale, device | 输入均为包含两个元素的列表，<br />第一元素代表ref，第二个元素代表prb                                                        |
    | 2 | x         | Nodes的三维坐标 | (2, bs*nn, pm)          | -                                                        |                                                                                                                             |
    | 3 | edges     | Edges输入向量   | (2, 2, bs*(n^2-3n+1))   | -                                                        | 列表中每一个元素<br />形如[tensor1, tensor2]这样的方式构建；<br />代表[rows, cols];<br />只记录邻接矩阵中的上三角矩阵中的边 |
    | 4 | edge_attr | 边属性向量      | (2, bs*(n^2-3n+1), 5)   | -                                                        |                                                                                                                             |
    | 5 | node_mask | 掩模padding节点 | (2, bs, nn)             | charges                                                  | 动态padding;<br />未使用；                                                                                                  |
    | 6 | edge_mask | 掩模不存在边    | (2, bs*(n^2-3n+1), 1)   |                                                          |                                                                                                                             |
    | 7 | n_nodes   | 节点数量        | (2, 1)                  |                                                          |                                                                                                                             |
    | 8 | label     | 标签            | (1)                     |                                                          | actives还是decoys                                                                                                           |

### 2.5.3 数据输入
1. 数据预处理
    使用molecule-3d-similarity/test.ipynb文件中的代码，对DUD-E的数据进行预处理，对每一个靶点得到具有`(name, smiles, label, charges, position, edges)`列的CSV文件。
2. DataUp
    `waiting...`
3. Dataset
    `waiting...`
4. DataLoader
    `waiting...`

# Question
1. sudo: add-apt-repository: command not found
    https://linuxconfig.org/sudo-apt-add-repository-command-not-found
    > sudo apt install software-properties-common

2. 
    