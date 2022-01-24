#export PATH=/Applications/MATLAB_R2022a.app/bin/:$PATH
alias matlab="/Applications/MATLAB_R2022a.app/bin/matlab -nodisplay -nojvm -r" # -nodesktop"

romsfile=$1
swanfile=$2
cwd=$(pwd)

echo "Converting ROMS matfile $romsfile to netCDF4...."

matlab "mat2roms_mw('${cwd}/${romsfile}', 'roms_grid.nc'); exit"

echo "Converting SWAN matfile $swanfile to netCDF4..."

matlab "mat2roms_mw('$swanfile', 'swan_grid.nc'); exit"

echo "Converting SWAN netCDF4 grid to .grid and .bot files"

matlab "roms2swan('swan_grid.nc'); exit"

fval=-0.000727
echo "Imposing f-plane on ROMS grid" 


