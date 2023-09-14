"""ML utilities."""

import torch
from torchmetrics import Accuracy
accuracy = Accuracy(task="binary")

def calculate_masked_accuracy(
    input: torch.Tensor,
    target: torch.Tensor,
    mask: torch.Tensor,
) -> torch.Tensor:
    """
    Calculate the accuracy of the model's predictions.

    Only considers the predictions for the masked residues, 
    i.e. where the mask is 1.0.

    Parameters
    ----------
    input : torch.Tensor
        The model's predictions.
    target : torch.Tensor
        The target values.
    mask : torch.Tensor
        The mask tensor.

    Returns
    -------
    torch.Tensor
        The accuracy.

    """
    # Get indexes where mask is 1.0
    mask_idx = torch.where(mask == 1)
    acc = accuracy(input[mask_idx], target[mask_idx])
    return acc
    


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
    def __init__(
        self, 
        subtype: str = "bce", # or "logits"
    ):
        self.subtype = subtype
        super(MaskedBinaryCrossEntropy, self).__init__()
        
    
    def forward(
        self,
        input, 
        target,
        mask,
    ):
        if self.subtype == "bce":
            
            loss = torch.nn.functional.binary_cross_entropy(
                input=torch.flatten(input),
                target=torch.flatten(target),
                weight=torch.flatten(mask),
                reduction="none",
            )

        elif self.subtype == "logits":
            # BCE with logits
            loss = torch.nn.functional.binary_cross_entropy_with_logits(
                input=torch.flatten(input),
                target=torch.flatten(target),
                weight=torch.flatten(mask),
                reduction="none",
            )

        result = torch.sum(loss) / int(torch.sum(mask))
        return result
