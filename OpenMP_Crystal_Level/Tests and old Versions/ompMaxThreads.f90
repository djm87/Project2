
program ompMaxThreads 
implicit none 
integer :: omp_get_max_threads 
print *, omp_get_max_threads() 
end program 
