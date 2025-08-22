wget https://github.com/sylabs/singularity/releases/download/v${SINGULARITY_VERSION}/singularity-ce-${SINGULARITY_VERSION}.tar.gz
tar -xzf singularity-ce-${SINGULARITY_VERSION}.tar.gz
cd singularity-ce-${SINGULARITY_VERSION}
./mconfig
make -j -C builddir
sudo make -j -C builddir install
