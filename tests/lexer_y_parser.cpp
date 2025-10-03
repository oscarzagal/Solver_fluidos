// Referencia: https://llvm.org/docs/tutorial/MyFirstLanguageFrontend/index.html

#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

/*-----------------------------------------------------------------------------
                                    Lexer
-----------------------------------------------------------------------------*/

enum Token {

    tok_eof        = -1,

    // Comandos
    tok_def        = -2,
    tok_extern     = -3,

    // Primario
    tok_identifier = -4,
    tok_number     = -5,

};

static std::string IdentifierStr; // Se llena si tok_identifier
static double NumVal;             // Se llena si tok_number

    static int LastChar = ' ';
    static int gettok() {

    // Saltar espacios en blanco
    while (isspace(LastChar)) {
        LastChar = getchar();
    }

    // Identificacion de [a-zA-Z][a-zA-Z0-9]* (tokens que inicien con letras)
    if (isalpha(LastChar)) {
        IdentifierStr = LastChar;
        while (isalnum(LastChar = getchar())) {
            IdentifierStr = LastChar;
        }

        if (IdentifierStr == "def") return tok_def;
        if (IdentifierStr == "extern") return tok_extern;

        return tok_identifier;
    }

    // Identificacion de [0-9.]+ (tokens que son numeros)
    if (isdigit(LastChar || LastChar == '.')) {
        std::string NumStr;
        do {
            NumStr += LastChar;
            LastChar = getchar();
        } while (isdigit(LastChar) || LastChar == '.');

        // Convierte de string a un valor numerico
        // La funcion "c_str()" supongo que convierte un std::string a un string
        // de C
        NumVal = strtod(NumStr.c_str(), 0);
        return tok_number;
    }

    // Saltar comentarios
    if (LastChar == '#') {

        do {
            LastChar = getchar();
        } while (LastChar != EOF && LastChar != '\n' && LastChar != '\r');

        if (LastChar != EOF) return gettok();
    }

    // Si la entrada no hace match con los casos de arriba o si se llega al
    // final del archivo
    if (LastChar == EOF) return tok_eof;

    // Retorno del caracter desconocido
    int ThisChar = LastChar;
    LastChar = getchar();
    return ThisChar;

}

int main() {

    gettok();

    return 0;
}

/*-----------------------------------------------------------------------------
                                Fin Lexer
  -----------------------------------------------------------------------------*/
