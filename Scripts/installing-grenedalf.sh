##installing grenedalf

# 1 Load the necessary software packages via the command


module load htslib/1.10.2  cmake/3.18.4  bzip2/1.0.8   autoconf/2.69   automake/1.16.3

# 2. Download grenedalf and start the installation

cd ~/software

git clone --recursive https://github.com/lczech/grenedalf.git
cd grenedalf
make clean
make

# After this installation the grenedalf executable was created under the grendedalf/bin directory:
cd bin
grenedalf
/grenedalf --help