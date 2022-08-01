###
 # @Author: haoqiang haoqiang@mindrank.ai
 # @Date: 2022-07-29 03:14:05
 # @LastEditors: haoqiang haoqiang@mindrank.ai
 # @LastEditTime: 2022-07-29 03:16:39
 # @FilePath: /work-home/molecule-3d-similarity/scripts/install-acpc.sh
 # @Description: Install ACPC
 # 
 # Copyright (c) 2022 by haoqiang haoqiang@mindrank.ai, All Rights Reserved. 
### 


sudo apt update
sudo apt install software-properties-common -y
sudo add-apt-repository ppa:avsm/ppa # may get message that 'ERROR: '~avsm' user or team does not exist', please try again
sudo apt update
sudo apt install gnuplot-x11 autoconf opam make gcc patch -y

rm -rf ~/.opam
opam init --disable-sandboxin # for LXC container

opam install ACPC


