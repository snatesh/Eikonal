#!/bin/bash

rm -rf $1
mkdir $1
mv $3 $1
mv $2 $1
mv $4 $1
lrztar -zf $1
stat --printf="%s\n" $1.tar.lrz
#rm -rf $1 $1.tar.lrz
