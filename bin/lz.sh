#!/bin/bash

mkdir bench
mv bench.coeffs bench/
mv bench.stl bench/
lrztar -zf bench/
val=$(stat --printf="%s" bench.tar.lrz)
rm -rf bench
echo $val
