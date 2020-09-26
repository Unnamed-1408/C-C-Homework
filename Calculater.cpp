#include <iostream>
#include <string>
#include<cstring>
#include<algorithm>
#include <sstream>
#include <math.h>



using namespace std;
void push(char a);
void pop();
void storage(string inital);
int priority(char in);
bool check_priority(char input , char yuan);
bool isclear();
void meet_brackets();
int find_comma(string a);
string pre_subtraction(string a , string b);
string pre_add(string a , string b);
int find_comma(string* a);
void add_back_zero(string* a);
void add_front_zero(string* a);
string start_add(string a , string b , int Aptr , int Bptr , string addoutput);
string addition(string a , string b);
string start_subtraction(string a , string b , int Aptr , int Bptr , string addoutput);
string substraction(string a , string b);
string ans_pop();
void calculate();
string split_flow(char a);
void ans_psuh(string a);
string martix(const char *first ,const char *second ,int lengtha ,int lengthb);
string multiply_addition(int **New , int lengtha , int lengthb);
int find_comma(string* a);
void remove_comma(string* a , int ptr);
void add_comma(string* a ,int ptr);
string ready_multiply(string a , string b);
string pre_multiply(string a , string b);
void insert_zero(string* a , int ptr);
void pre_process_input_string(string* a);
string start_divide(string a , string b);
bool isvariables(string a);
void storage_variables(string a);
string start(string a);
int search_variables(string a);
void replace_variables(string*a , int variables_ptr_temp , int last , int* ptr);
void pre_process_variables_replace(string* a);
void restore();
void re_input(string* a);
int ptr_function_right(string a , int left);
int search_function(string a);
string pow(string a);
string select_function(int name , string a);
void function_replace(string* a , int left , int right , string to_replace);
void check_brackets(string a);
bool has_function(string a);
void function_start(string* a);
bool find_letter(string a);
string sqrt(string a);
void function_storage();
void check_input(string a);
void trim_output(string* a);
void help();
string sin(string a);
string cos(string a);
string tan(string a);
string abs(string a);
void trim_blank(string* a);

int flag_brackets = -1;
string output = "";
char stack[100];
int stackptr = 0;//zhi xiang yi ge kong zhan zheng
string anstack[100];
int ans_stackptr = 0;
int error_flag = -1;

string variables[100];
string variables_values[100];
int variables_ptr = 0;

string call_function_name[10];
const int function_num = 6;

int main()
{
    function_storage();
    help();
    while(1){
    string bc = "";
    cout <<"[in]";
    getline(cin,bc);
/*    pre_process_input_string(&bc);
    storage(bc);
    cout << output;
    calculate();
    cout<<endl;
    cout<<anstack[ans_stackptr-1];
    */
    if(bc=="cl" || bc=="clear"){
        system("clear");
    }
    else{

        if(bc=="h" || bc=="help"){
            help();
    }
    else{
        if(bc=="q" || bc=="quit"){
            cout << "bye" << endl;
            exit(0);
        }
        else{
            trim_blank(&bc);
            string temp = start(bc);
            if(temp != ""){
                trim_output(&temp);
                cout << "==> " << temp << endl;
            }
            else{
            }
        }
    }
    }
    }
}

void trim_blank(string* a){
    int ptr = 0;
    string temp = "";
    while(ptr != (*a).length()){
        if(a->at(ptr) != ' '){
            temp += a->at(ptr);
        }
        ptr++;
    }
    (*a) = temp;
}


void function_storage(){
    call_function_name[0] = "pow";
    call_function_name[1] = "sqrt";
    call_function_name[2] = "sin";
    call_function_name[3] = "cos";
    call_function_name[4] = "tan";
    call_function_name[5] = "abs";
}

void restore(){
    flag_brackets = -1;
    output = "";
    stackptr = 0;
    ans_stackptr = 0;
    error_flag = -1;
}

string start(string a){
    if(isvariables(a)){

        int ptr = 0;
        while(ptr != a.length()){
            if(a.at(ptr) == '='){
                break;
            }
            ptr++;
        }

        if(search_variables(a.substr(0,ptr)) == -1){
            storage_variables(a);
        }
        else{
            if((start(a.substr(ptr+1))) != ""){

         //   variables[variables_ptr] = a.substr(0,ptr);
        //    cout << variables[variables_ptr];
            if((start(a.substr(ptr+1))).at(0) == '-'){
                variables_values[search_variables(a.substr(0,ptr))] =    start(a.substr(ptr+1)) + ')';
                variables_values[search_variables(a.substr(0,ptr))] = '0'+variables_values[search_variables(a.substr(0,ptr))];
                variables_values[search_variables(a.substr(0,ptr))] = '('+variables_values[search_variables(a.substr(0,ptr))];
                cout << "variable " << variables[search_variables(a.substr(0,ptr))] << " has been updated!"<<endl;
        //        cout<<variables_values[variables_ptr];
            }
            else{
                variables_values[search_variables(a.substr(0,ptr))] = '(' + start(a.substr(ptr+1)) + ')';
                cout << "variable " << variables[search_variables(a.substr(0,ptr))] << " has been updated!"<<endl;
            }
        //    cout<<variables_values[variables_ptr];
        //    cout<<variables_values[variables_ptr];
            variables_ptr++;
            }
            else{

            }

/*            string b = a.substr(ptr+1);
            pre_process_variables_replace(&b);

            b = start(b);
            check_input(b) ;

            if(error_flag == -1){
                if(b.at(0) == '-'){
                    variables_values[search_variables(a.substr(0,ptr))] =    start(a.substr(ptr+1)) + ')';
                    variables_values[search_variables(a.substr(0,ptr))] = '0'+variables_values[search_variables(a.substr(0,ptr))];
                    variables_values[search_variables(a.substr(0,ptr))] = '('+variables_values[search_variables(a.substr(0,ptr))];
//                    cout << variables_values[search_variables(a.substr(0,ptr))];
                    cout << "variable " << variables[search_variables(a.substr(0,ptr))] << " has been updated!"<<endl;
                }
                else{
                    variables_values[search_variables(a.substr(0,ptr))] = '(' +  b + ')';
//                	cout << variables_values[search_variables(a.substr(0,ptr))];
                    cout << "variable " << variables[search_variables(a.substr(0,ptr))] << " has been updated!"<<endl;
                }
            }
            else{
                re_input(&b);
            }*/
        //    cout << variables_values[search_variables(a.substr(0,ptr))];
        //    cout<<variables_values[variables_ptr];
        //    variables_ptr++;
        }
    }
    else{
        pre_process_input_string(&a);
        check_input(a);

        if(error_flag == -1){
            storage(a);
//        	cout << output;
            calculate();
//        	cout<<endl;
            restore();
            return anstack[0];
        }
        else{
            re_input(&a);
            restore();
            error_flag = -1;
            return "";
        }

    }
    return "";
}

bool isvariables(string a){
    int ptr = 0;
    while(ptr != a.length()){
        if(a.at(ptr) == '='){
            return true;
        }
        ptr++;
    }
    return false;
}

void storage_variables(string a){
    int ptr = 0;
    while(ptr != a.length()){
        if(a.at(ptr) == '='){
            break;
        }
        ptr++;
    }
    if((start(a.substr(ptr+1))) != ""){

    variables[variables_ptr] = a.substr(0,ptr);
//    cout << variables[variables_ptr];
    if((start(a.substr(ptr+1))).at(0) == '-'){
        variables_values[variables_ptr] =    start(a.substr(ptr+1)) + ')';
        variables_values[variables_ptr] = '0'+variables_values[variables_ptr];
        variables_values[variables_ptr] = '('+variables_values[variables_ptr];
//        cout<<variables_values[variables_ptr];
    }
    else{
    variables_values[variables_ptr] = '(' + start(a.substr(ptr+1)) + ')';
    }
//    cout<<variables_values[variables_ptr];
//    cout<<variables_values[variables_ptr];
    variables_ptr++;
    }
    else{

    }
}

void pre_process_variables_replace(string* a){
    int ptr = 0;
    int last = 0;
//    int length = (*a).length();

    while(ptr < (*a).length()){
//        cout << *a;
        if(((*a).at(ptr) >='a' && (*a).at(ptr) <='z') || ((*a).at(ptr) >='A' && (*a).at(ptr) <='Z')){
            last = ptr;
            while((*a).at(ptr) != '*' &&(*a).at(ptr) != '/'&&(*a).at(ptr) != '+'&&(*a).at(ptr) != '-'&&(*a).at(ptr) != ')' && (*a).at(ptr)!='('){
                ptr++;
                if(ptr == (*a).length()){
                    break;
                }
            }
//            cout << (*a);
            if(ptr < (*a).length() && (*a).at(ptr) == '('){
                return;
            }

            int temptr = search_variables((*a).substr(last , ptr-last));
            if(temptr == -1){
                error_flag = 1 ;
                return;
            }
            else{
                replace_variables(a , temptr , last , &ptr);
            }
        }
        ptr++;
    }
}

bool find_letter(string a){
    int ptr = 0;
    while(ptr != a.length()){
        if(isalpha(a.at(ptr))){
            return true;
        }
        ptr++;
    }
    return false;
}

void re_input(string* a){
    cout << "What you input has something wrong!"<<endl;
    cout << "The reason maybe no variables has found"<<endl;
    cout << "Plz check your input and input again"<<endl;

}

void replace_variables(string*a , int variables_ptr_temp , int last , int* ptr){
    (*a) = (*a).substr(0 , last) + variables_values[variables_ptr_temp] + (*a).substr(*ptr);
    (*ptr) += last+(variables_values[variables_ptr_temp]).length()-1;
}


int search_variables(string a){
    int ptr = 0;
    while(ptr != variables_ptr){
        if(variables[ptr] != a){
            ptr++;
        }
        else{
            return ptr;
        }
    }
    return -1;
}

void pre_process_input_string(string* a){
    int ptr = 0;
    while(ptr != (*a).length()){

        if((*a).at(ptr) == '-'){

            if(ptr == 0){
                insert_zero(a , 0);
            }

            else if((*a).at(ptr-1) == '('){
                insert_zero(a,ptr);
            }
        }
        ptr++;
    }
//    cout << *a;
    check_brackets(*a);
    if(error_flag != 1){
//        cout << (*a);
        while(find_letter(*a)){

            if(error_flag == 1){
                break;
            }
//        pre_process_variables_replace(a);
        if(has_function(*a)){
            function_start(a);
    }
    else{
    }

        if(find_letter(*a) == false){
            return;
        }
//        cout << *a;
        pre_process_variables_replace(a);

        }
    }
}

void check_input(string a){
    int ptr = 0;
    while(ptr != a.length()){
        if(a.at(ptr) == '*' || a.at(ptr) == '/' || a.at(ptr) == '+' || a.at(ptr) == '-' || a.at(ptr) == '=' ){
            while(ptr != a.length()-1 && a.at(ptr+1) == ' '){
                ptr++;
            }
            if(ptr !=a.length()-1 && ((a.at(ptr+1)) == '*' || (a.at(ptr+1) == '/' || (a.at(ptr+1)) == '+' || (a.at(ptr+1)) =='-' || (a.at(ptr+1)) =='='))){
                error_flag = 1;
//                re_input(&a);
                return;
            }
        }
        ptr++;
    }
    ptr = 0;
    while(ptr != a.length()){
        if(a.at(ptr) == '*' || a.at(ptr) == '/' || a.at(ptr) == '+' || a.at(ptr) == '-' || a.at(ptr) == '=' ){
            while(ptr != a.length()-1 && a.at(ptr+1) == ' '){
                ptr++;
            }
            if(ptr ==a.length()-1 || a.at(ptr+1) == ')'){
                error_flag = 1;
//                re_input(&a);
                return;
            }
        }
        ptr++;
    }
    ptr = 0;
    while(ptr != a.length()){
        if(isalnum(a.at(ptr)) || a.at(ptr) == '(' || a.at(ptr) == ')' || a.at(ptr) == '+' || a.at(ptr) == '-' || a.at(ptr) == '*' || a.at(ptr) == '/' || a.at(ptr) == '.'){
            ptr++;
        }
        else{
            error_flag = 1;
//            re_input(&a);
            return;
        }
    }
    return;
}

void check_brackets(string a){
    int ptr = 0;
    int count = 0;
    while(ptr != a.length()){
        if(a.at(ptr) == '('){
            count++;
        }
        if(a.at(ptr) == ')'){
            count--;
        }
        ptr++;
    }
    if(count != 0){
        error_flag = 1;
//        re_input(&a);
    }
    else{
        return;
    }
}


void insert_zero(string* a , int ptr){

    (*a) = (*a).substr(0 , ptr) + '0' + (*a).substr(ptr);

}

/*void pre_prosess(string a){
    int ptr = 0;
    while(ptr != a.length()){

        if(a.at(ptr) == '-'){
            a.insert(ptr-1,"0");
        }

        int temp = ptr-1;
        while(a.at(temp) == ' '){
            if(isdigit(a.at(temp))){


            }
        }
        ptr++;
    }
}

*/


void storage(string inital){
    int ptr = 0;
    while(ptr != inital.length()){

        if(isdigit(inital.at(ptr)) || inital.at(ptr) == '.' || inital.at(ptr) == ' '){
            output+=inital.at(ptr);
        }
        else{
            output+=" ";
            push(inital.at(ptr));
        }
        ptr++;
    }
    output+=" ";
    while(isclear() == false){
        pop();
    }

}

void push(char a){
        if(a == ')'){
            while(flag_brackets == -1){
                meet_brackets();
            }
            flag_brackets = -1;
        }
        else{
              while(check_priority(a , stack[stackptr-1]) == false && stack[stackptr-1] != '('){
                 pop();
             }
             stack[stackptr] = a;
             stackptr++;
        }
}

void pop(){

    if(stack[stackptr-1] != '('){
        output+=stack[stackptr-1];
        output+=" ";
    }
    else {
        flag_brackets = 1;
    }
    stackptr--;
}


bool check_priority(char input , char yuan){
    if(priority(input) > priority(yuan)){
        return true;
    }
    else {
        return false;
    }
}

int priority(char in){
    switch (in) {
        case '+' : return 0; break;
        case '-' : return 0; break;
        case '*' : return 1; break;
        case '/' : return 1; break;
        case '(' : return 2; break;
        default: return -1;
    }
}

bool isclear(){
    if(stackptr != 0){
        return false;
    }
    else {
        return true;
    }

}

void meet_brackets(){
    pop();
}

void calculate(){

    int ptr = 0;
    int last = 0;
    while(ptr != output.length()){

        while(output.at(ptr) == ' '){
            ptr++;
            if(ptr == output.length()){
                return;
            }
        }
        last = ptr;

        if( isdigit(output.at(ptr)) ){

            while(isdigit(output.at(ptr)) || output.at(ptr) == '.'){
                  ptr++;
                  if(ptr == output.length()){
                      break;
                  }
            }

            anstack[ans_stackptr] = output.substr(last , ptr-last);
//            cout << anstack[ans_stackptr] <<" ";
            ans_stackptr++;
        }
        else{
            ans_psuh(split_flow(output.at(ptr)));
            ptr++;

        }

    }
}

string split_flow(char a){
    switch (a) {
        case '+' : return pre_add(ans_pop() , ans_pop());
        case '-' : return pre_subtraction(ans_pop() , ans_pop());
        case '*' : return pre_multiply(ans_pop() , ans_pop());
        case '/' : return start_divide(ans_pop() , ans_pop());
    }

}

string ans_pop(){

    ans_stackptr--;

    return anstack[ans_stackptr];

}

void ans_psuh(string a){
    anstack[ans_stackptr] = a;
    ans_stackptr++;
}

string pre_add(string a , string b){
    if(a.at(0) == '-' && b.at(0) == '-'){
        string temp = "";
        temp = addition(a.substr(1) , b.substr(1));
        return '-'+temp;
    }
    if(a.at(0) == '-' || b.at(0) == '-'){

        if(a.at(0) == '-'){
            return pre_subtraction(b , a.substr(1));
        }
        else {
            return pre_subtraction(a , b.substr(1));
        }

    }
    else {
        return addition(a , b);
    }

}


string pre_subtraction(string a , string b){

    if(a.at(0) == '-' && b.at(0) == '-'){
        string temp = "";
        temp = substraction(a.substr(1) , b.substr(1));
        return '-'+temp;
    }
    if(a.at(0) == '-' || b.at(0) == '-'){

        if(a.at(0) == '-'){
            return pre_add(a , '-'+b);
        }
        else {
            return pre_add(a , b.substr(1));
        }
    }
    else {
        return substraction(a , b);
    }

}

string pre_multiply(string a , string b){

    if((a.at(0) == '-' && b.at(0) == '-') || (a.at(0) != '-' && b.at(0) != '-')){
        if(a.at(0) == '-' && b.at(0) == '-'){
            return ready_multiply(a.substr(1) , b.substr(1));
        }
        else {
            return ready_multiply( a , b );
        }
    }
    else {

            if(a.at(0) == '-'){
                return '-' + ready_multiply(a.substr(1) , b);
            }
            else {
                return '-' + ready_multiply(a , b.substr(1));
            }
    }


}

string ready_multiply(string a , string b){

    int Aptr = find_comma(&a);
    int Bptr = find_comma(&b);

    int o_lengtha = a.length();
    int o_lengthb = b.length();

    remove_comma(&a , Aptr);
    remove_comma(&b , Bptr);
    int last = o_lengtha-Aptr + o_lengthb-Bptr ;

    int lengtha = a.length();
    int lengthb = b.length();

    reverse(b.begin(),b.end());

    const char *first = a.c_str();
    const char *second = b.c_str();
//    int **New = martix(first , second , lengtha , lengthb);
    string output_multiply = martix(first , second , lengtha , lengthb);
    add_comma(&output_multiply , last);
    if(output_multiply.at(0) == '0'){
        output_multiply = output_multiply.substr(1);
    }

    return output_multiply;

}



//加法计算器
string addition(string a , string b){
    int Aptr = find_comma(&a);
    int Bptr = find_comma(&b);
    string addoutput = "";

    return start_add(a , b , Aptr , Bptr , addoutput);
}


int find_comma(string* a){

    int temptr = 0;

    while(temptr != (*a).length()){

        if((*a).at(temptr) == '.'){
            return temptr;
        }
        temptr++;
    }

    (*a) += '.';
    return temptr;

}

string start_add(string a , string b , int Aptr , int Bptr , string addoutput){
//    int tempa = Aptr;
//    int tempb = Bptr;
    int forward = 0;

    if(a.length()-Aptr > b.length()-Bptr){
        while(a.length()-Aptr != b.length()-Bptr){
            add_back_zero(&b);
        }
    }
    else {
        while(a.length()-Aptr != b.length()-Bptr){
            add_back_zero(&a);
        }
    }//back zero add

    if(Aptr > Bptr){
        while(Aptr != Bptr){
            add_front_zero(&b);
                Bptr++;
        }
    }
    else {
        while(Aptr != Bptr){
            add_front_zero(&a);
            Aptr++;
        }
    }

//    if(a.length() == b.length()){
//        cout << "well done!";
//    }

    int temp = a.length()-1;

    while(temp != -1){

        if(a.at(temp) == '.'){
            addoutput = '.' + addoutput;
            temp--;
        }
        else{
         int abc = a.at(temp)-'0' + b.at(temp) - '0' + forward;
         int now = abc%10;
         addoutput = (char)(now + '0') + addoutput;
         forward = abc/10;
         temp--;
    }


}
    if(forward != 0){
        addoutput = (char)(forward + '0')+addoutput;
    }

    return addoutput;
}


void add_back_zero(string* a){
    (*a)+='0';
}

void add_front_zero(string* a){
    (*a) = '0' + *a;
}


//end



//减法计算器
string substraction(string a , string b){
    int Aptr = find_comma(&a);
    int Bptr = find_comma(&b);
    string addoutput = "";

    return start_subtraction(a , b , Aptr , Bptr , addoutput);
}


/*int find_comma(string* a){

    int temptr = 0;

    while(temptr != (*a).length()){

        if((*a).at(temptr) == '.'){
            return temptr;
        }
        temptr++;
    }

    (*a) += '.';
    return temptr;

}
*/
string start_subtraction(string a , string b , int Aptr , int Bptr , string addoutput){
//    int tempa = Aptr;
//    int tempb = Bptr;
    int forward = 0;

    if(a.length()-Aptr > b.length()-Bptr){
        while(a.length()-Aptr != b.length()-Bptr){
            add_back_zero(&b);
        }
    }
    else {
        while(a.length()-Aptr != b.length()-Bptr){
            add_back_zero(&a);
        }
    }//back zero add

    if(Aptr > Bptr){
        while(Aptr != Bptr){
            add_front_zero(&b);
                Bptr++;
        }
    }
    else {
        while(Aptr != Bptr){
            add_front_zero(&a);
            Aptr++;
        }
    }

/*    if(a.length() == b.length()){
        cout << "well done!";
    }
*/
    int temp = a.length()-1;

    while(temp != -1){

        if(a.at(temp) == '.'){
            addoutput = '.' + addoutput;
            temp--;
        }
        else{
         int abc = a.at(temp)-'0' - b.at(temp) + '0' + forward;
         if(abc < 0){
             addoutput = (char)(10 + abc + '0') + addoutput;
             forward = -1;
             temp--;
         }
         else{
            forward = 0;
            addoutput = (char)(abc + '0') + addoutput;
             temp--;
        }
    }

}
    if(forward != 0){
        if(forward < 0){
            return '-' + substraction(b,a);
        }
        else{
        addoutput = (char)(forward + '0')+addoutput;
        }
    }

    return addoutput;
}

//end



//乘法计算器
string martix(const char *first ,const char *second ,int lengtha ,int lengthb){
    int *martix[lengthb] ;
    int count = 0 ;

    while(count != lengthb){

        int temp = lengtha;
        int *tempm = new int[lengtha+1] ();
        int forward = 0;
        while(temp != -1){

            if(temp != 0){

            tempm[temp] = ((first[temp-1]-'0')*(second[count]-'0') + forward+0) % 10;

            forward = (((first[temp-1]-'0')*(second[count]-'0') + forward) - tempm[temp])/10;
        }
            else {
                tempm[temp] = forward ;

            }
       //     cout << tempm[temp] ;

            temp--;

        }
        forward = 0;
        martix[count] = tempm;
        count++;
    }

    return multiply_addition( martix , lengtha , lengthb);
}

string multiply_addition(int *New[] , int lengtha , int lengthb){
    int forward = 0;
    string temp = "";
    int movex = lengtha;
    int movey = 0;

    while(movex != 0 || movey != lengthb){
        int tempx = movex;
        int tempy = movey;
        int total = 0;
        int tempa = 0;
        while(tempx != lengtha+1 && tempy != lengthb){

            total += New[tempy][tempx];
            tempx++;
            tempy++;

        }
        if(movex != 0){
           movex--;
        }
        else{
            movey++;

        }
        tempa = total;
        total = total+forward;
        forward = ((tempa +forward)-(tempa +forward)%10)/10;
        total = total%10;
        temp+=(total+'0');

    }

    reverse(temp.begin(),temp.end());
    return temp;

}

/*int find_comma(string* a){

    int temptr = 0;

    while(temptr != (*a).length()){

        if((*a).at(temptr) == '.'){
            return temptr;
        }
        temptr++;
    }

    (*a) += '.';
    return temptr;

}
*/
void remove_comma(string* a , int ptr){
    *a = (*a).substr(0,ptr)+(*a).substr(ptr+1);
}

void add_comma(string* a ,int ptr){
    *a = (*a).substr(0,(*a).length()-ptr+2) + '.' +(*a).substr((*a).length()-ptr+2);
}


//end


//除法计算器
string start_divide(string a , string b){
    double A = atof(a.c_str());
    double B = atof(b.c_str());

    double divide_ans = A/B;
    ostringstream os;
    os << std::fixed <<divide_ans;
    return os.str();
}
//end



//函数
bool has_function(string a){
    int ptr = 0 ;
    int last = 0;

    while(ptr < (a).length()){
        if(((a).at(ptr) >='a' && (a).at(ptr) <='z') || ((a).at(ptr) >='A' && (a).at(ptr) <='Z')){
            last = ptr;
            while((a).at(ptr) != '*' &&(a).at(ptr) != '/'&&(a).at(ptr) != '+'&&(a).at(ptr) != '-' && (a).at(ptr) != ')' && a.at(ptr) != '('){
                ptr++;
                if(ptr == (a).length()){
                    break;
                }
            }
            if(ptr < a.length() && a.at(ptr) == '('){
//                cout << ptr_function_right(a , ptr);
                return true;
            }

        }
        ptr++;
    }
    return false;
}

void function_start(string* a){
    int ptr = 0 ;
    int last = 0;
    string name = "";
    string parameter = "";
    while(ptr < (*a).length()){
        if(((*a).at(ptr) >='a' && (*a).at(ptr) <='z') || ((*a).at(ptr) >='A' && (*a).at(ptr) <='Z')){
            last = ptr;
            while((*a).at(ptr) != '*' &&(*a).at(ptr) != '/'&&(*a).at(ptr) != '+'&&(*a).at(ptr) != '-' && (*a).at(ptr) != ')' && (*a).at(ptr) != '('){
                ptr++;
                if(ptr == (*a).length()){
                    break;
                }
            }
            if(ptr < (*a).length() && (*a).at(ptr) == '('){
//                cout << ptr_function_right(a , ptr);
                parameter  = (*a).substr(ptr+1 , ptr_function_right(*a,ptr)-ptr-1);
                name = (*a).substr(last , ptr-last);

                if(search_function(name) != -1){
                    function_replace(a , last , ptr_function_right(*a,ptr) ,select_function(search_function(name) , parameter));
                }
                else{
                    error_flag = 1;
                }
            }

        }
        ptr++;
    }

}


void function_replace(string* a , int left , int right , string to_replace){
    (*a) = (*a).substr(0,left) + to_replace + (*a).substr(right+1);
}

string select_function(int name , string a){
    switch (name) {
        case  0 : return pow(a);
        case 1 : return sqrt(a);
        case 2: return sin(a);
        case 3 : return cos(a);
        case 4: return tan(a);
        case 5: return abs(a);
    }

}

string pow(string a){
//    cout << a;
    bool temp = false;
    int temptr = 0;
    while(temptr != a.length()){
        if(a.at(temptr) == ','){
            temp = true;
        }
        temptr++;
    }
    if(a == "" || temp == false){
        cout << "error from pow , the output/variables will be set as 0" << endl;
        return "0";
    }

    int ptr = 0;
    int last = 0;
    int count = 0;
    while(ptr != a.length()){

        if(a.at(ptr) == '('){
            count++;
        }
        if(a.at(ptr) == ')'){
            count--;
        }

        if(count == 0 && a.at(ptr) == ','){
            break;
        }
        ptr++;
   }


    if(ptr == a.length()){
//        cout << "pow invaild!"<<endl;
//        cout << "To protect.....exiting..."<<endl;
//        exit(0);
        error_flag = 1;
    }
    else{
        string num = start(a.substr(0,ptr));
        string count = start(a.substr(ptr+1));
        string temp = to_string(pow(atof(num.c_str()) , atof(count.c_str())));
//        cout << temp;
        if(temp == "inf"){
            cout << "Sorry that , the number is too large... it's beyond double...pow is crash" << endl;
            error_flag = 1;
            return "";
        }
        if(temp.at(0) == '-'){
            return "(0-"+temp.substr(1)+')';
        }
        else{
            return '('+temp+')';
        }
    }
}

string sqrt(string a){
    string temp = to_string(sqrt(atof(start(a).c_str())));
    return temp;
}

/*
string pow_process(string a , string b){
    string output_pow = "";
    while((int)atof(b.c_str()) != 0){

    }
}
*/

string sin(string a){

        double tempa = (atof(start(a).c_str()))/180;
        while(tempa >=2){
            tempa-=2;
        }

        string temp = to_string(sin(tempa*M_PI));
        return temp;

}

string cos(string a){

        double tempa = (atof(start(a).c_str()))/180;
        while(tempa >=2){
            tempa-=2;
        }

        string temp = to_string(cos(tempa*M_PI));
        return temp;

}

string tan(string a){

        double tempa = (atof(start(a).c_str()))/180;
        while(tempa >=2){
            tempa-=2;
        }

        string temp = to_string(tan(tempa*M_PI));
        return temp;

}


string abs(string a){

    if(start(a).at(0) == '-'){
        return start(a).substr(1);
    }
    else{
        return a;
    }

}

int ptr_function_right(string a , int left){
    int ptr = left;
    int count = 0;
    while(ptr != a.length()){

        if(a.at(ptr) == '('){
            count++;
        }
        if(a.at(ptr) == ')'){
            count--;
        }

        if(count == 0){
            return ptr;
        }

        ptr++;

    }
    return 0;
}

int search_function(string a){
    int ptr = 0;
    while(ptr != function_num){
        if(call_function_name[ptr] != a){
            ptr++;
        }
        else{
            return ptr;
        }
    }
    return -1;
}

void trim_output(string* a){
//  cout << *a;
    bool temp = false;
    int temptr = 0;
    while(temptr != (*a).length()){
        if((*a).at(temptr) == '.'){
            temp = true;
        }
        temptr++;
    }



    int ptr = (*a).length()-1;
    if(temp == true){
    while(ptr != -1){
       if((*a).at(ptr) == '0' ){
            ptr--;
        }
        else{
            break;
        }
    }

    if(ptr != -1 && (*a).at(ptr) == '.'){
        ptr--;
    }

    (*a) = (*a).substr(0,ptr+1);
    }

    ptr = 0;
    int flag = 0;

    while(ptr != (*a).length()){

        if((*a).at(ptr) == '-'){
            flag = 1;
            ptr++;
        }

        if((*a).at(ptr) == '0'){
            ptr++;
        }
        else{
            break;
        }

    }
    if(flag == 1){
        (*a) = '-' + (*a).substr(ptr);
    }
    else {
        (*a) = (*a).substr(ptr);
    }
    if((*a).length() == 0 ||((*a).length()==1 && flag == 1)){
        (*a) = "0";
    }
    ptr = 0;
    while (ptr != (*a).length()) {
        if((*a).at(ptr) == '.'){
            break;
        }
        ptr++;
    }

    if(ptr -1 < 0 || isdigit((*a).at(ptr-1)) == false ){
        if(ptr-1>=0){
            (*a) = (*a).substr(0,ptr)+'0'+(*a).substr(ptr);
        }
        else{
            (*a) = '0'+(*a).substr(0);
        }
    }
}


void help(){
    cout << "*****************************************高精度浮点计算器************************************" << endl;
    cout << "                                      author:Xie yuezhen                                   " << endl;
    cout << "                    此计算器暂时只支持 + - * / 四则基本运算，以及一些简单函数的调用                 " << endl;
    cout << "Sample:" << endl;
    cout << "1+2     2*3    1*(2+3)   sqrt(2)    pow(2,4)    pow(pow(2,4),4)" << endl;
    cout<<endl;
    cout << "目前支持的函数有:" << endl;
    cout << "pow  sqrt abs sin cos tan   除了abs支持高精度计算" << endl;
    cout << "此计算器支持函数的叠用" << endl;
    cout << "sqrt(sqrt(sqrt(10)))  sqrt(x)"<<endl;
    cout << "但是函数的调用均不支持高精度,均为调用math库进行运算,有上限double,以及涉及到除法及函数的时候会有少许精度损失" << endl;
    cout<<endl;
    cout << "此计算器支持变量的存储" << endl;
    cout << "Sample:"<<endl;
    cout << "x=1 y=1 xyzjjoq=1"<<endl;
    cout << "注意：变量不能包含数字及其他符号，只能由字母组成" <<endl;
    cout << "支持变量的变化" << endl;
    cout << "x=x+1 y=x+1 x=sqrt(4)" << endl;
    cout << "******************************-introduction end-**********************************************"<<endl;
    cout<<endl;
    cout << "**********************************warning*****************************************************"<<endl;
    cout <<"                                    非法输入                                                     "<<endl;
    cout << "1+-2   1+   2+(   x+=2                                           "<<endl;
    cout << "在遇到1+-2的时候，可以将其换成1+(-2)输入"<<endl;
    cout <<"*************************************-end-******************************************************"<<endl;
    cout << R"(    可输入"h"或者"help"再次调出,可输入"q"或者"quit"退出,在linux环境下可用"cl"或者"clear"进行清屏 )"<<endl;
    cout << R"(                           h/help/q/quit/cl/clear均作为保留字使用                         )"<<endl;
    cout << "************************************************************************************************"<<endl;
}





























