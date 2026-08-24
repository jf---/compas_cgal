#pragma once

#include "audit_depletion_witness_2.h"

AuditDepletionWitness2 apply_audit_segment_depletion_to_trial(
    const Stock2& authority,
    Stock2& trial,
    const AuditSegmentMotion2& motion,
    const AuditPolicy2& policy);
AuditDepletionWitness2 apply_audit_circle_depletion_to_trial(
    const Stock2& authority,
    Stock2& trial,
    const AuditCircleMotion2& motion,
    const AuditPolicy2& policy);
AuditDepletionWitness2 apply_audit_arc_depletion_to_trial(
    const Stock2& authority,
    Stock2& trial,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy);
AuditDepletionWitness2 apply_audit_plunge_depletion_to_trial(
    const Stock2& authority,
    Stock2& trial,
    const AuditVerticalPlunge2& motion,
    const AuditPolicy2& policy);

bool validate_audit_segment_depletion_witness(
    const Stock2& authority,
    const Stock2& post,
    const AuditSegmentMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDepletionWitness2& witness);
bool validate_audit_circle_depletion_witness(
    const Stock2& authority,
    const Stock2& post,
    const AuditCircleMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDepletionWitness2& witness);
bool validate_audit_arc_depletion_witness(
    const Stock2& authority,
    const Stock2& post,
    const AuditArcMotion2& motion,
    const AuditPolicy2& policy,
    const AuditDepletionWitness2& witness);
bool validate_audit_plunge_depletion_witness(
    const Stock2& authority,
    const Stock2& post,
    const AuditVerticalPlunge2& motion,
    const AuditPolicy2& policy,
    const AuditDepletionWitness2& witness);
