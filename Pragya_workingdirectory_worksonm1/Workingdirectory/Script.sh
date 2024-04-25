for i in `seq 1 8`
do
    mkdir -p ASSYMETRIC/COLLISION/SEED${i}

    echo "#!/bin/bash -l" >> Asym_Coll${i}.sbat
    echo "" >>Asym_Coll${i}.sbat

    echo "#SBATCH --job-name=Asym${i}" >>Asym_Coll${i}.sbat
    echo "#SBATCH --account=lappi" >>Asym_Coll${i}.sbat
    echo "#SBATCH --partition=small" >>Asym_Coll${i}.sbat 
    echo "#SBATCH --time=16:00:00" >>Asym_Coll${i}.sbat
    echo "#SBATCH --nodes=1" >>Asym_Coll${i}.sbat
    echo "#SBATCH --cpus-per-task=40" >>Asym_Coll${i}.sbat

    echo "" >>Asym_Coll${i}.sbat

    echo "module add fftw/3.3.10-mpi-omp gsl" >>Asym_Coll${i}.sbat

    echo "" >>Asym_Coll${i}.sbat

    echo "./ASym_SUNcSolver.exe -Nx 128 -Ny 128 -Nz 1024 -SEED ${i} -o ASSYMETRIC/COLLISION/SEED${i}/" >>Asym_Coll${i}.sbat
    
    sbatch Asym_Coll${i}.sbat

done





for i in `seq 1 8`
do

    mkdir -p ASSYMETRIC/LEFT_NUCLEUS/SEED${i}

    echo "#!/bin/bash -l" >> Asym_Left${i}.sbat
    echo "" >>Asym_Left${i}.sbat

    echo "#SBATCH --job-name=Asym_Left${i}" >>Asym_Left${i}.sbat
    echo "#SBATCH --account=lappi" >>Asym_Left${i}.sbat
    echo "#SBATCH --partition=small" >>Asym_Left${i}.sbat 
    echo "#SBATCH --time=16:00:00" >>Asym_Left${i}.sbat
    echo "#SBATCH --nodes=1" >>Asym_Left${i}.sbat
    echo "#SBATCH --cpus-per-task=40" >>Asym_Left${i}.sbat

    echo "" >>Asym_Left${i}.sbat

    echo "module add fftw/3.3.10-mpi-omp gsl" >>Asym_Left${i}.sbat

    echo "" >>Asym_Left${i}.sbat

    echo "./ASymLeft_SUNcSolver.exe -Nx 128 -Ny 128 -Nz 1024 -SEED ${i} -o ASSYMETRIC/LEFT_NUCLEUS/SEED${i}/" >>Asym_Left${i}.sbat

    sbatch Asym_Left${i}.sbat
done





for i in `seq 1 8`
do

    mkdir -p ASSYMETRIC/RIGHT_NUCLEUS/SEED${i}

    echo "#!/bin/bash -l" >> Asym_Right${i}.sbat
    echo "" >>Asym_Right${i}.sbat

    echo "#SBATCH --job-name=Asym_Right${i}" >>Asym_Right${i}.sbat
    echo "#SBATCH --account=lappi" >>Asym_Right${i}.sbat
    echo "#SBATCH --partition=small" >>Asym_Right${i}.sbat 
    echo "#SBATCH --time=16:00:00" >>Asym_Right${i}.sbat
    echo "#SBATCH --nodes=1" >>Asym_Right${i}.sbat
    echo "#SBATCH --cpus-per-task=40" >>Asym_Right${i}.sbat

    echo "" >>Asym_Right${i}.sbat

    echo "module add fftw/3.3.10-mpi-omp gsl" >>Asym_Right${i}.sbat

    echo "" >>Asym_Right${i}.sbat

    echo "./ASymRight_SUNcSolver.exe -Nx 128 -Ny 128 -Nz 1024 -SEED ${i} -o ASSYMETRIC/RIGHT_NUCLEUS/SEED${i}/" >>Asym_Right${i}.sbat
    
    sbatch Asym_Right${i}.sbat
done





for i in `seq 1 8`
do
    mkdir -p SYMMETRIC/COLLISION/SEED${i}

    echo "#!/bin/bash -l" >> Coll${i}.sbat
    echo "" >>Coll${i}.sbat

    echo "#SBATCH --job-name=Sym${i}" >>Coll${i}.sbat
    echo "#SBATCH --account=lappi" >>Coll${i}.sbat
    echo "#SBATCH --partition=small" >>Coll${i}.sbat 
    echo "#SBATCH --time=16:00:00" >>Coll${i}.sbat
    echo "#SBATCH --nodes=1" >>Coll${i}.sbat
    echo "#SBATCH --cpus-per-task=40" >>Coll${i}.sbat

    echo "" >>Coll${i}.sbat

    echo "module add fftw/3.3.10-mpi-omp gsl" >>Coll${i}.sbat

    echo "" >>Coll${i}.sbat

    echo "./SUNcSolver.exe -Nx 128 -Ny 128 -Nz 1024 -SEED ${i} -o SYMMETRIC/COLLISION/SEED${i}/" >>Coll${i}.sbat
    
    sbatch Coll${i}.sbat
done





for i in `seq 1 8`
do

    mkdir -p SYMMETRIC/LEFT_NUCLEUS/SEED${i}

    echo "#!/bin/bash -l" >> Left${i}.sbat
    echo "" >>Left${i}.sbat

    echo "#SBATCH --job-name=Left${i}" >>Left${i}.sbat
    echo "#SBATCH --account=lappi" >>Left${i}.sbat
    echo "#SBATCH --partition=small" >>Left${i}.sbat 
    echo "#SBATCH --time=16:00:00" >>Left${i}.sbat
    echo "#SBATCH --nodes=1" >>Left${i}.sbat
    echo "#SBATCH --cpus-per-task=40" >>Left${i}.sbat

    echo "" >>Left${i}.sbat

    echo "module add fftw/3.3.10-mpi-omp gsl" >>Left${i}.sbat

    echo "" >>Left${i}.sbat

    echo "./LeftSUNcSolver.exe -Nx 128 -Ny 128 -Nz 1024 -SEED ${i} -o SYMMETRIC/LEFT_NUCLEUS/SEED${i}/" >>Left${i}.sbat

    sbatch Left${i}.sbat 
done





for i in `seq 1 8`
do

    mkdir -p SYMMETRIC/RIGHT_NUCLEUS/SEED${i}

    echo "#!/bin/bash -l" >> Right${i}.sbat
    echo "" >>Right${i}.sbat

    echo "#SBATCH --job-name=Right${i}" >>Right${i}.sbat
    echo "#SBATCH --account=lappi" >>Right${i}.sbat
    echo "#SBATCH --partition=small" >>Right${i}.sbat 
    echo "#SBATCH --time=16:00:00" >>Right${i}.sbat
    echo "#SBATCH --nodes=1" >>Right${i}.sbat
    echo "#SBATCH --cpus-per-task=40" >>Right${i}.sbat

    echo "" >>Right${i}.sbat

    echo "module add fftw/3.3.10-mpi-omp gsl" >>Right${i}.sbat

    echo "" >>Right${i}.sbat

    echo "./RightSUNcSolver.exe -Nx 128 -Ny 128 -Nz 1024 -SEED ${i} -o SYMMETRIC/RIGHT_NUCLEUS/SEED${i}/" >>Right${i}.sbat
    
    sbatch Right${i}.sbat
done
