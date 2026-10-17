void BuildSimpleRMS(SimParameters *, SimBoxRange ,  size_t,  DoDeFunc *, DoDeInfo *, MPI_Comm );
void ReconstructRMS(SimParameters *, SimBoxRange ,  size_t,  DoDeFunc *, DoDeInfo *, MPI_Comm );
void refinempirms(void **, ptrdiff_t *, ptrdiff_t, DoDeFunc *, DoDeInfo *, MPI_Comm, GridInfo *);
void BroadCastDDFuncAndDDInfo(SimParameters *);
