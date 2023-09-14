

import torch


class MaskedMSELoss(torch.nn.Module):
    def __init__(self):
        super(MaskedMSELoss, self).__init__()
    
    def forward(
        self,
        input, 
        target,
        mask,
    ):
        diff2 = (torch.flatten(input) - torch.flatten(target)) ** 2.0 * torch.flatten(mask)
        result = torch.sum(diff2) / torch.sum(mask)
        return result


class MaskedBinaryCrossEntropy(torch.nn.Module):
    def __init__(self):
        super(MaskedBinaryCrossEntropy, self).__init__()
    
    def forward(
        self,
        input, 
        target,
        mask,
    ):
        # BCE with logits
        loss = torch.nn.functional.binary_cross_entropy_with_logits(
            input=torch.flatten(input),
            target=torch.flatten(target),
            weight=torch.flatten(mask),
            reduction="none",
        )
        result = torch.sum(loss) / int(torch.sum(mask))
        return result