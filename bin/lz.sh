#!/bin/bash

mkdir -p $1
mv $3 $1
mv $2 $1
lrztar -zf $1
stat --printf="%s\n" $1.tar.lrz
rm -rf $1 $1.tar.lrz
