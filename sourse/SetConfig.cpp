#include "./SetConfig.h"

using namespace std;

string ChToRem = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz _;:,=\\/%АБВГДЕЁЖЗИЙКЛМНОПРСТУФХЦЧШЩЪЫЬЭЮЯабвгдеёжзийклмнопрстуфхцчшщъыьэюя";

void RemoveChars(string &source)
{
   string result = "";
   bool inside = false;
   for (unsigned int i = 0; i < source.length(); i++)
   {
      if (source[i] == '\'')
         inside = !inside;

      if (inside && source[i] != '\'')
      {
         result += source[i];
      }
   }
   source = result;
}

void RemoveChars(string &source, const string &chars)
{
   string result = "";
   for (unsigned int i = 0; i < source.length(); i++)
   {
      bool foundany = false;
      for (unsigned int j = 0; j < chars.length() && !foundany; j++)
      {
         foundany = (source[i] == chars[j]);
      }
      if (!foundany)
      {
         result += source[i];
      }
   }
   source = result;
}

void GetParam(string &St, int &Param)
{
   cout << St << endl;
   RemoveChars(St, ChToRem);
   Param = stoi(St);
};

void GetParam(string &St, double &Param)
{
   cout << St << endl;
   RemoveChars(St, ChToRem);
   Param = stod(St);
};

void GetParam(string &St, string &Param)
{
   cout << St << endl;
   RemoveChars(St);
   Param = St;
};

string d_to_string(double &St)
{
   string str = to_string(St);
   str.erase(str.find_last_not_of('0') + 1, std::string::npos);
   return str;
}
