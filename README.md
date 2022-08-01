<!--
 * @Author: haoqiang haoqiang@mindrank.ai
 * @Date: 2022-07-19 03:41:19
 * @LastEditors: haoqiang haoqiang@mindrank.ai
 * @LastEditTime: 2022-07-29 03:22:30
 * @FilePath: /work-home/molecule-3d-similarity/README.md
 * @Description: 
 * 
 * Copyright (c) 2022 by haoqiang haoqiang@mindrank.ai, All Rights Reserved. 
-->
# 1 Overview
| 序号 | 名称 | 来源 | 相关文献 | 安装步骤 |  
| -- | -- | -- | -- | -- |
| 1 | ShapeProtrudeDist | rdkit | - | - |
| 2 | ShapeTanimotoDist | rdkit | - | - |
| 3 | ShapeTverskyIndex | rdkit | - | - |
| 4 | SC score | rdkit | @Yang2020 | - |
| 5 | Gobbi_Pharm2D | rdkit | - | - |
| 6 | USR | oddt | @Ballester2007 | - |
| 7 | USRCAT | oddt | @Schreyer2012 | - |
| 8 | Electroshape | oddt | @Armstrong2010 | - |
| 9 | ACPC | [riken](http://www2.riken.jp/zhangiru/software.html) | @Berenger2014 | [Installation](#ACPC) |
|  |  |  |  | - |
|  |  |  |  | - |
|  |  |  |  | - |
|  |  |  |  | - |
|  |  |  |  | - |
|  |  |  |  | - |
|  |  |  |  | - |
|  |  |  |  | - |
|  |  |  |  | - |
|  |  |  |  | - |
|  |  |  |  | - |


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


# Question
1. sudo: add-apt-repository: command not found
    https://linuxconfig.org/sudo-apt-add-repository-command-not-found
    > sudo apt install software-properties-common

2. 
    