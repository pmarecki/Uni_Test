#include <cstdio>
#include <cstdlib>

#define KMAX 10000
#define THREADS 32
#define BLOCKS 32
#define LoopMAX 2000
#define klucz 137

//czas GPU: user  0m0.336s
//czas CPU: user  0m44.523s


//Funkcje uruchamialne na GPU z CPU ---> "__global__"
__global__ void funkcja(int *miejsce) {
  int indeks = threadIdx.x + blockIdx.x * blockDim.x;
//  printf("th=%i block=%i dane_GPU[%i]=%i\n",
//      threadIdx.x, blockIdx.x, indeks, miejsce[indeks]);
  int liczba = miejsce[indeks];

  for(int i=0; i<LoopMAX; ++i)
    for(int j=0; j<LoopMAX; ++j)
      liczba = (liczba + i * j) % klucz;

  miejsce[indeks] = liczba;
}


void funkcja_CPU(int *miejsce) {
  for(int indeks = 0; indeks < THREADS * BLOCKS; ++indeks) {
    int liczba = miejsce[indeks];
    for(int i=0; i<LoopMAX; ++i)
      for(int j=0; j<LoopMAX; ++j)
        liczba = (liczba + i * j) % klucz;
    miejsce[indeks] = liczba;
  }

}



int main(void) {
  int threads_per_block = THREADS;
  int blocks = BLOCKS;

  int *dane;
  dane = (int*) malloc(KMAX * 4);     //alokacja CPU

  int *dane_GPU;
  cudaMalloc(&dane_GPU, KMAX * 4);    //alokacja GPU

  //wygenerowanie liczb (CPU)
  for(int i=0; i<KMAX; ++i)
    dane[i] = i;

  //kopiowanie danych CPU --> GPU
  //syntax: (cel, zrodlo, ilosc byteow, flaga)
  //przesylanie powrotne (dane, dane_GPU, KMAX * 4, cudaMemcpyDeviceToHost)
//  cudaMemcpy(dane_GPU, dane, KMAX * 4, cudaMemcpyHostToDevice);
//  funkcja<<<threads_per_block, blocks>>>(dane_GPU);
//  cudaDeviceSynchronize();              //oczekiwanie na koniec obliczen GPU


  funkcja_CPU(dane);


  free(dane);  //dealokacja pamieci (CPU)
  cudaFree(dane_GPU);
  return 0; 
}
