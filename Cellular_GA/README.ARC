Kódovanie: ISO-8859-2

Informácie k testovaciemu programu testOneMax
---------------------------------------------

Je pripravený program testOneMax, ktorý implementuje rie¹enie OneMax problému so
zvolenou då¾kou chromozómu.

Preklad programu:

make

Spustenie programu a nápoveda:

testOneMax -h

Då¾ku chromozómu je mo¾né nastavi» v zdrojovom programe testOneMax.cc,
nastavením premennej chrLength na po¾adovanú hodnotu.

Popis vypracovania projektu pre predmet ARC/EVO
-----------------------------------------------

1. Úvod

Podµa ¹peciálneho zadania bol implementovaný celulárny model genetického
algoritmu, ktorý pracuje paralelne. Bolo pou¾ité rozhranie MPI, implementované v
MPICH.


2. Mapovanie na procesory

Pre mapovanie na procesory bolo pou¾ité ¹achovnicové mapovanie. Ka¾dý procesor
vlastní maticu jedincov. Napríklad pre 9 procesorov dostaneme nasledujúce
rozdelenie:

         +--------------+
         | P0 | P1 | P2 |
         +----+----+----+
         | P3 | P4 | P5 |
         +----+----+----+
         | P6 | P7 | P8 |
         +----+----+----+

  Obr. 1: Mapovanie na procesory

Poèet procesorov je obmedzený na mocninu 2, aby bolo mo¾né vytvori» rozumné
mapovanie. Veµkos» matice jedincov v genetickom algoritme v¹ak nie je obmedzený.
Ka¾dý procesor takto vlastní èas» globálnej mrie¾ky jedincov.


3. Komunikácia

V celulárnom modele ka¾dý jedinec komunikuje so svojimi susedmi z pevne
definovaného okolia. V paralelnej verzii je toto okolie obmedzené na vzdialenos»
jedna (manhattan), t.j. jedinec si vyberá pre krí¾enie v¾dy svojich susedov
sprava, zµava, zhora, alebo zdola. Prebiehajú teda tieto typy komunikácie:
 
 1. Vysladnie k ¾iadosti o chromozóm k susedovi
 2. Príjem ¾iadosti o chromozóm
 3. Reakcia na ¾iados» o chromozóm od suseda -- preposlanie chromozómu
 4. Preposlanie výsledkov evolúcie na lokálnej mrie¾ke hlavnému uzlu (rank 0),
    po skonèení jednej generácie
 5. Preposlanie správy z hlavného uzlu v¹etkým ostatným procesorom s
    informáciou o ïal¹om pokraèovaní evolúcie. Sú tri mo¾nosti: pokraèuj,
    koniec, po¹li mi najlep¹ieho jedinca.

Správy 1. -- 3. prebiehajú medzi susednými uzlami a teda medzi susediacimi
procesormi. Správy 4. a 5. prebiehajú medzi hlavným procesorom a v¹etkými
ostatnými procesormi.

Pre komunikáciu boli pou¾ité dve zásady:
 
 1. Príjem správ je v¾dy neblokujúci
 2. Vysielanie správ je blokujúce,

a to z toho dôvodu, ¾e ak by obidve varianty boli blokujúce mohlo by dochádza» k
deadlocku. Napríklad v prípade podµa obr. 1 je mo¾né, ¾e procesory P1, P2, P4, a
P5 vytvoria cyklus ¾iadostí o jedinca: P1 --> P2 --> P5 --> P4 --> P1. Ak na
procesoroch nebol pred vyslaním ¾iadosti o chromozóm u¾ dopredu spustený
neblokujúci príjem tejto ¾iadosti, do¹lo by k deadlocku.

Dôsledkom tohoto systému sú dve veci:

 1. Ka¾dý procesor musí poèas vykonávania genetického algoritmu sledova» príjem
    ¾iadostí (príjem je zaèatý pri ¹tarte lokálnej verzie) od svojich
    susedov a v prípade po¾iadavky na chromozóm tento chromozóm vysla» k
    susedovi.
 2. Procesor, ktorý po¾aduje nejaký chromozóm od svojho suseda, musí pred
    vyslaním po¾iadavky zaèa» neblokujúci príjem chromozómu a následne v
    aktívnej smyèke testova» vybavenie ¾iadosti. V aktívnej smyèke ïalej musí
    procesor testova» aj ¾iadosti od ostatných procesorov.

Týmto spôsobom sa zabráni deadlocku pri komunikácii medzi procesormi. Nevýhodou
je v¹ak, ¾e pri spustení paralelnej verzie na jednoprocesorovom systéme
dochádza k nadmernému za»a¾eniu procesoru pri aktívnom èakaní, èo brzdí ostatné
procesy vo vykonávaní programu. Takto je mo¾né v celom systém pozorova»
spomalenie oproti sekvenènej verzii.

Pri spustení na viacerých procesoroch v¹ak dochádza k zrýchleniu a pri testovaní
som ani raz nepozoroval spomalenie.


4. Namerané výsledky

V archíve je pripravený skript ARCtest.sh, ktorý spustí niekoµkokrát genetický
algoritmus podµa zadaných parametrov. Výsledkom je tabuµka obsahujúca poèet
procesorov, èas výpoètu, zrýchlenie a efektivitu výpoètu. Pred spustením skriptu
je nutne správne nastavi» práva a parametre skriptu mpi.

Algoritmus som testoval na dvoch variantách, s pou¾itím jednoduchého problému
OneMax:

 1. 
   chromosomeLength=350 (då¾ka chromozómu v bitoch)
   xsize=300 ysize=300 gen=100 NUM_PROCS="1 4 9 16 25 36" ./ARCtest.sh

   #Processors     Time [s]     Speedup         Effectivity
   1               815.156      1.00000         1.00000
   4               254.035      3.20883          .80220
   9               129.181      6.31018          .70113
   16              183.881      4.43306          .27706
   25              119.365      6.82910          .27316
   36              65.7563      12.39662         .34435

 2.
   chromosomeLength=500
   xsize=600 ysize=600 gen=50 NUM_PROCS="1 4 9 16 25 36" ./ARCtest.sh

   #Processors/cores       Time [s]        Speedup         Effectivity
   1                       2320.6         1.00000         1.00000
   4                       652.631        3.55576         .88894
   9                       366.246        6.33617         .70401
   16                      329.66         7.03937         .43996
   25                      184.991        12.54439        .50177
   36                      156.4          14.83759        .41215


Je vidie», ¾e s rastúcim poètom procesorov klesá efektivita paralelného
spracovania k nule. Predpokladám, ¾e je to spôsobené komunikáciou medzi
procesormi a aktívnym èakaním pri neblokujúcom spracovaní správ.
