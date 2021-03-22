#! /usr/bin/env bash

ls Data/scr/ | awk -F "_" '{if($2>y){y=$2}  } ;END{print y}'
