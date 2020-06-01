#ifndef BEAST_DECODERID_H
#define BEAST_DECODERID_H

#include <string>

enum class DecoderID {
  ERROR,
  BEAST,
  BSD,
  KTKL,
  ORDERED_STATISTICS,
  VITERBI,
};

static std::string idToString(DecoderID id) {
  switch (id) {
    case DecoderID::BEAST:return "beast";
    case DecoderID::BSD:return "bsd";
    case DecoderID::ORDERED_STATISTICS:return "ordered-statistics";
    case DecoderID::KTKL:return "ktkl";
    case DecoderID::VITERBI:return "viterbi";
    default:return "";
  }
}
static DecoderID stringToId(const std::string &name) {
  if (name == "beast") {
    return DecoderID::BEAST;
  } else if (name == "bsd") {
    return DecoderID::BSD;
  } else if (name == "ktkl") {
    return DecoderID::KTKL;
  } else if (name == "ordered-statistics") {
    return DecoderID::ORDERED_STATISTICS;
  } else if (name == "viterbi") {
    return DecoderID::VITERBI;
  } else {
    return DecoderID::ERROR;
  }
}
#endif //BEAST_DECODERID_H
