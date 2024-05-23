
void temp(){
	std::vector<int> arr {1, 2, 3, 4, 5};
	const int alen = 5;
	for (int i = 0; i < alen; i++){
		cout << i << ": ";
		int jstart;

		switch(i){
		case 0:
			jstart = 0;
			break;
		case alen - 1:
			jstart = alen + 1;
			break;
		default:
			jstart = i + 1;
		}

		cout << "jstart: " << jstart << endl;
		for (int j = jstart; j <= i + 1; j++){
			cout << arr[j] << " ";
		}
		cout << endl << "*";
	}
}