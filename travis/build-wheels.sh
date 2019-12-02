#!/bin/bash
set -e -x

shopt -s expand_aliases
# pick the oldest version of python to install cmake:
# PYDIR=`ls -d /opt/python/cp3*/bin | head -1`
/opt/python/cp34-cp34m/bin/pip install --upgrade cmake
cp /io/travis/fake_cmake.sh /bin/cmake && chmod +x /bin/cmake
#alias cmake="/opt/python/cp34-cp34m/bin/python -c 'import cmake; cmake.cmake()'"
## install up-to-date pybind11 headers
# Download jq and place it where we can find it:
curl -s -L https://github.com/stedolan/jq/releases/download/jq-1.6/jq-linux64 -o /bin/jq
chmod +x /bin/jq
# Grab an up-to-date pybind11:
pybind11_url=`curl -s "https://api.github.com/repos/pybind/pybind11/tags"|jq -r '.[0].tarball_url'`
curl -s -L "${pybind11_url}" -o pybind11.tar.gz
mkdir -p pybind11/build
tar -xzf pybind11.tar.gz -C pybind11 --strip-components 1
cd pybind11/build && cmake .. -DPYBIND11_TEST=OFF && make install


# install dependencies
for PYBIN in /opt/python/cp3*/bin; do
	"${PYBIN}/pip" install --upgrade -r /io/travis/requirements.txt
done
# perform the actual building
for PYBIN in /opt/python/cp3*/bin; do
	"${PYBIN}/pip" wheel /io/ -w wheelhouse/
done

# bundle external shared libraries into the wheels
for wheel in wheelhouse/*.whl; do
	auditwheel repair "$wheel" --plat $PLAT -w /io/wheelhouse/
done

# Install packages and test
for PYBIN in /opt/python/cp3*/bin/; do
	"${PYBIN}/pip" install brille --no-index -f /io/wheelhouse
#	(cd "$HOME"; ${PYBIN}/nosetests" brille)
done
