#! /usr/bin/env bash

scp -i ./andres2.pem aws* ubuntu@$1:~

ssh -i ./andres2.pem ubuntu@$1


