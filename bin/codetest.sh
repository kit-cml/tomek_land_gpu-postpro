cd ..
make clean all
cd bin
./drug_sim -input_deck input_deck_example.txt -hill_file control_me/control.csv -init_file control_me/init_cell1.csv
