#ifndef BEAST_INCLUDE_EXCEPTIONS_H_
#define BEAST_INCLUDE_EXCEPTIONS_H_
#include <exception>

class InvalidFileException : std::exception {
  [[nodiscard]] const char* what() const noexcept override {
    return "Invalid file!";
  }
};
#endif //BEAST_INCLUDE_EXCEPTIONS_H_
