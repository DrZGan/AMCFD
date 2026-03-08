find -name "*.o" -exec rm -f '{}' \;
find -name "*.f90~" -exec rm -f '{}' \;
find -name "*.mod" -exec rm -f '{}' \;
find -name "*.log" -exec rm -f '{}' \;
find -name "cluster_main*" -exec rm -f '{}' \;

cd result
find -mindepth 1 -maxdepth 1 -type d -exec rm -rf '{}' \;