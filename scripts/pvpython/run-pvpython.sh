#!/bin/bash

# $1 = paraview (pv)python macro file template
# $2 = basename of XDMF file(s)

function print_usage {
  echo "`basename $0`"
  echo "Usage: `basename $0` pvpyton_macro.template adam_solution_file_basename"
  echo
  echo "Example:"
  echo "  `basename $0` z-slice-extract.template sphere-shock-000009954"
  echo
}

if [ $# -eq 0 ] ; then
  print_usage
  exit
fi

template=`echo $1 | sed -e 's/\.template//'`
for xdmf in $2*.xdmf ; do
   it=`echo $xdmf | awk -F '-' '{print $3}' | sed -e 's/^0*//' | sed -e 's/\.xdmf//'`
   if [ -f "z-slice-$it.vtm" ]; then
      echo "$it already processed"
   else
      echo $xdmf
      cat $template.template | sed -e "s/AAAIT/-$it/" | sed -e "s/BBFILENAMEXDMF/$xdmf/" > $template-$it.py
      pvpython --force-offscreen-rendering $template-$it.py
   fi
done
