#!/bin/bash

cd ../build && make install
cd ../python && python3 eikonal_solve.py
