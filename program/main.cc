/*
 * Copyright (C) 2016- Satoko Yamamoto, Natsuki Hosono, Yoko Funato, Junichiro Makino
 *
 * This code is licensed under MIT license found in the LICENSE file..
 * 
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <math.h>
#include <iostream>
#include <getopt.h>
#include <vector>

typedef std::vector<double> doublevector;

#include "class.h"
#include "sampledisk.h"
#include "userdefineddisk.h"
#include "wendland-62.h"
#include "rings.h"
#include "disks.h"

int main(int argc, char* argv[]){
	return DiskCreator<Factory> (Input (argc, argv)).make_disk();
}


