#!/bin/bash

dir=$1
for file in "$dir"/*
do
  sed -i -E "s/^(>NODE_[0-9]+)_.+$/\1/" "$file"
done