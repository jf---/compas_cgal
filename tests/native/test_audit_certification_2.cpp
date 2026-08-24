void audit_certification_contract_gate();
void audit_certification_verdicts_gate();
void audit_certification_refinement_gate();
void audit_certification_falsifiers_gate();
void audit_certification_identity_gate();

void audit_certification_gate()
{
    audit_certification_contract_gate();
    audit_certification_verdicts_gate();
    audit_certification_refinement_gate();
    audit_certification_falsifiers_gate();
    audit_certification_identity_gate();
}
