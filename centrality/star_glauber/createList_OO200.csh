#!/bin/csh -f
# Make LIST/tree.<type>.list from the O+O Glauber trees (spherical), as createList.csh does for Ru+Ru
set types = ( "default" "small" "large" "smallXsec" "largeXsec" "gauss" "smallNpp" "largeNpp" )

foreach type ($types)
  set treelist = "./LIST/tree.$type.list"
  echo "Make $treelist..."
  ls -1 ./output/fastglaubermc_OO_200GeV_${type}_spherical_*.root > $treelist
end

# smallTotal/largeTotal/lowrw/highrw re-use the default trees (as in createList.csh)
cd LIST
ln -vsf tree.default.list tree.smallTotal.list
ln -vsf tree.default.list tree.largeTotal.list
ln -vsf tree.default.list tree.lowrw.list
ln -vsf tree.default.list tree.highrw.list
ls -lrta
cd -
