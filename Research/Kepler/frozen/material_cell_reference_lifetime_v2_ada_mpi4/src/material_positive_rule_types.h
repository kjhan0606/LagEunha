#ifndef MATERIAL_POSITIVE_RULE_TYPES_H
#define MATERIAL_POSITIVE_RULE_TYPES_H
// Native read-only rule payload. Mode0 evaluates the original full rule.
struct PositiveRule {
    bool loaded=false,active=false;
    double anchor_time=0;
    uint64_t anchor_state=0,generation=0,global_points=0;
    std::string directory,binding;
    std::vector<uint64_t>offset,count,row_offset;
    std::vector<uint32_t>index;
    std::vector<double>weight;
    std::vector<R>residual,rounding;
    std::vector<std::array<R,16>>diagnostic;
    std::vector<uint8_t>compressed;
};
#endif
