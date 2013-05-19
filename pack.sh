#!/bin/bash
#
#   pack.sh
#
#   Pack all the files for modeldb/senselab
#
#       Copyright (C) 2013  Lukas Solanka <l.solanka@sms.ed.ac.uk>
#       
#       This program is free software: you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation, either version 3 of the License, or
#       (at your option) any later version.
#       
#       This program is distributed in the hope that it will be useful,
#       but WITHOUT ANY WARRANTY; without even the implied warranty of
#       MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#       GNU General Public License for more details.
#       
#       You should have received a copy of the GNU General Public License
#       along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
pack_dir=GridCellModel

rm -f grid_cell_model/*.pyc
mkdir -p grid_cell_model/output
rm -f grid_cell_model/output/*

cd ..
tar -c -z --exclude .git* -f $pack_dir.tar.gz $pack_dir
