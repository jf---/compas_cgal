#pragma once

#include "partition_certificate.h"

#include <string>
#include <vector>

std::string encode_string_sequence(
    const std::vector<std::string>& values);

std::vector<std::string> decode_string_sequence(
    const std::string& encoded);

void finalize_event_partition(
    EventPartitionCertificate2& certificate);

VerifiedEventPartition2 verify_event_partition(
    const EventPartitionCertificate2& certificate);

EventPartitionCertificate2 mutate_certificate_record(
    const EventPartitionCertificate2& certificate,
    const std::string& mutation);
